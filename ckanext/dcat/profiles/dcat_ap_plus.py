import json
import os
from decimal import Decimal, DecimalException
import requests
from rdflib import term, URIRef, BNode, Literal, Graph
import ckantoolkit as toolkit

# from ckan.lib.munge import munge_tag
import logging

# NOTE: We import dataclasses from 'dcat_4c_ap' (local copy) instead of the
# official pip packages because the server runs Python 3.7, while the
# official packages require Python >3.9. Once the server is upgraded,
# revert to: from dcat_ap_plus.datamodel.dcat_ap_plus import ...
from ckanext.dcat.profiles.dcat_4c_ap import (Agent,
                                              Concept,
                                              Dataset,
                                              DataGeneratingActivity,
                                              DefinedTerm,
                                              Document,
                                              EvaluatedEntity,
                                              Entity,
                                              Identifier,
                                              LinguisticSystem,
                                              LegalResource,
                                              Standard,
                                              QualitativeAttribute,
                                              QuantitativeAttribute)
from .euro_dcat_ap import EuropeanDCATAPProfile
log = logging.getLogger(__name__)

from ckanext.dcat.utils import (
    resource_uri,
    DCAT_EXPOSE_SUBCATALOGS,
    DCAT_CLEAN_TAGS,
    publisher_uri_organization_fallback,
)
from .base import RDFProfile, URIRefOrLiteral, CleanedURIRef
from .base import (
    RDF,
    XSD,
    SKOS,
    RDFS,
    DCAT,
    DCT,
    ADMS,
    VCARD,
    FOAF,
    SCHEMA,
    NFDI,
    CHEMINF,  # this
    CHMO,  # this
    OBI,
    IAO,
    PROV,
    CHEBI,
    NMR,
    QUDT,
    NCIT,
    FIX,
    namespaces
)
from linkml_runtime.dumpers import RDFLibDumper
from linkml_runtime.utils.schemaview import SchemaView

# Module-level cache for SchemaView instances
# Key: schema_name, Value: SchemaView object
_SCHEMA_VIEW_CACHE = {}


class Helpers(object):
    """
    Mixin class containing shared helper methods for DCAT-AP+ and ChemDCAT-AP profiles.
    Handles data extraction, normalization, schema loading, and LinkML object instantiation.
    """

    # Class-level cache for PubChem CID lookups (persists across requests in a worker)
    _pubchem_cache = {}
    _CACHE_MAX_SIZE = 500

    def _get_schema_view(self, schema_name, local_filename, purl):
        """
        UNIVERSAL Schema Loader.
        Loads a SchemaView with a robust fallback strategy:
        1. Memory Cache
        2. Local File (Sibling 'schemas' folder)
        3. Remote PURL

        Args:
            schema_name: Unique key for caching (e.g., 'dcat_ap_plus')
            local_filename: Name of the YAML file (e.g., 'dcat_ap_plus.yaml')
            purl: Remote URL fallback

        Returns: SchemaView instance or None
        """
        # 1. Check Memory Cache
        if schema_name in _SCHEMA_VIEW_CACHE:
            return _SCHEMA_VIEW_CACHE[schema_name]

        schema_content = None
        source = ""

        # 2. Dynamic Path Calculation (Works for both profiles)
        # Assumes: profiles/profile.py and schemas/file.yaml are siblings
        current_dir = os.path.dirname(os.path.abspath(__file__))
        local_yaml_path = os.path.normpath(os.path.join(current_dir, "..", "schemas", local_filename))

        # 3. Try Local File
        try:
            if os.path.exists(local_yaml_path):
                sv = SchemaView(local_yaml_path, merge_imports=True)
                _SCHEMA_VIEW_CACHE[schema_name] = sv
                log.info(f"Schema '{schema_name}' loaded from local file and cached: {local_yaml_path}")
                return sv
        except Exception as e:
            log.error(f"Failed to parse local schema '{schema_name}' at {local_yaml_path}: {e}")

        # 4. Try Remote PURL
        if not schema_content:
            try:
                log.debug(f"Fetching schema '{schema_name}' from PURL: {purl}")
                resp = requests.get(purl, headers={"Accept": "application/yaml, text/yaml"}, timeout=10)
                resp.raise_for_status()
                schema_content = resp.text
                source = "remote PURL"
            except Exception as e:
                log.error(f"Failed to fetch schema '{schema_name}' from remote: {e}")
                return None

        # 5. Parse and Cache
        if schema_content:
            try:
                sv = SchemaView(schema_content, merge_imports=True)
                _SCHEMA_VIEW_CACHE[schema_name] = sv
                log.info(f"Schema '{schema_name}' loaded from {source} and cached.")
                return sv
            except Exception as e:
                log.error(f"Failed to parse schema '{schema_name}': {e}")
                return None

        return None

    def _get_pubchem_cid(self, inchi_key=None, smiles=None):
        """
        Fetches CID from PubChem with a class-level cache.
        """
        key = None
        if inchi_key:
            key = inchi_key.strip().upper()
        elif smiles:
            key = smiles.strip()

        if not key:
            return None

        if key in self._pubchem_cache:
            return self._pubchem_cache[key]

        cid = None
        try:
            suffix = f"inchikey/{key}" if inchi_key else f"smiles/{key}"
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{suffix}/cids/TXT"
            resp = requests.get(url, headers={"Accept": "text/plain"}, timeout=3)
            if resp.status_code == 200 and resp.text.strip():
                result = resp.text.strip().split("\n")[0]
                if result.isdigit():
                    cid = result
        except Exception:
            pass

        self._pubchem_cache[key] = cid

        # Simple rotation if cache gets too big
        if len(self._pubchem_cache) > self._CACHE_MAX_SIZE:
            keys_to_remove = list(self._pubchem_cache.keys())[:100]
            for k in keys_to_remove:
                del self._pubchem_cache[k]

        return cid


    def _get_authors(self, dataset_dict):
        """
        Parses the author string into a list of Agent objects.
        Handles various formats: "Last, First", "Last, Initial", or lists of names.
        TODO:
            * Search Service should normalize the authors/creators to an object with first name, last name and PID.
            * DCAT-AP plus needs better handling for this, see: https://github.com/nfdi-de/dcat-ap-plus/issues/84

        Args:
            dataset_dict: The dataset metadata dictionary.

        Returns:
            A list of Agent objects representing the creators.
        """
        creators = []
        author_string = dataset_dict.get("author")

        if not author_string or not isinstance(author_string, str):
            return creators

        fragments = [f.strip() for f in author_string.split(",") if f.strip()]
        if not fragments:
            return creators

        is_single_author = (
                len(fragments) == 2 and
                len(fragments[1]) > 2 and
                fragments[1].replace(".", "").isalpha() and
                " " not in fragments[0]
        )

        full_names = []
        if is_single_author:
            full_names.append(f"{fragments[0]} {fragments[1]}")
        else:
            for fragment in fragments:
                clean_frag = fragment.strip()
                if not clean_frag:
                    continue
                is_likely_initial = len(clean_frag) <= 2 and clean_frag.replace(".", "").isalpha()
                if is_likely_initial and full_names:
                    last_name = full_names.pop()
                    full_names.append(f"{last_name} {clean_frag}")
                else:
                    full_names.append(clean_frag)

        for name in full_names:
            if name.endswith("."):
                parts = name.split()
                if parts and len(parts[-1]) > 3:
                    name = name[:-1]
            name = " ".join(name.split())
            if name:
                creators.append(Agent(
                    name=name,
                    type=Concept(preferred_label='person', description='A human being.')
                ))

        return creators

    def _get_dataset_id(self, dataset_dict):
        """Constructs the canonical Dataset IRI."""
        if dataset_dict.get("doi"):
            return f"https://doi.org/{dataset_dict.get('doi')}"
        raw_id = dataset_dict.get("id", "").strip()
        return f"https://search.nfdi4chem.de/dataset/{raw_id}"

    def _get_other_ids(self, dataset_dict):
        """Constructs a list of Identifier objects."""
        raw_id = dataset_dict.get("id", "").strip()
        other_ids = [Identifier(notation=f"https://search.nfdi4chem.de/dataset/{raw_id}",
                                title="Search Service ID",
                                description="The id of this dataset within the NFDI4Chem Search Service "
                                            "(https://search.nfdi4chem.de/)")]
        if dataset_dict.get("doi"):
            other_ids.append(Identifier(notation=dataset_dict.get("doi"),
                                        title="DOI",
                                        description="The DOI of this dataset"))
        return other_ids

    def _get_compound_id(self, dataset_dict, dataset_id):
        """Resolves the Compound IRI (PubChem CID or local fragment)."""
        cid = self._get_pubchem_cid(dataset_dict.get("inchi_key"), dataset_dict.get("smiles"))
        if cid:
            return f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
        return f"{dataset_id}#sample_compound"

    def _get_description(self, dataset_dict):
        """Extracts and cleans the description."""
        desc = dataset_dict.get('notes')
        return desc.strip() if desc else 'No description'

    def _get_language(self, dataset_dict):
        """Normalizes language code and returns a LinguisticSystem object."""
        raw = (dataset_dict.get('language') or 'en').strip().lower()
        code = 'de' if raw in ('deutsch', 'german', 'de') else 'en'
        return LinguisticSystem(title=code, description=f"http://id.loc.gov/vocabulary/iso639-1/{code}")

    def _get_publisher(self, dataset_dict):
        """Extracts organization info and creates a Publisher Agent."""
        org = dataset_dict.get("organization") or {}
        org_name = org.get("title") or org.get("display_name") or org.get("name") or "Unknown Organization"
        return Agent(
            name=org_name,
            type=Concept(
                preferred_label='Academia/Scientific organisation',
                description='http://purl.org/adms/publishertype/Academia-ScientificOrganisation'
            )
        )

    def _get_license(self, dataset_dict, dataset_id):
        """Extracts license information and returns a list of LegalResource objects."""
        if not dataset_dict.get('license_title'):
            return []
        title = dataset_dict['license_title']
        license_id = dataset_dict.get('license_id')
        license_url = dataset_dict.get('license_url')
        url = f"{dataset_id}#license_notspecified" if (license_id == 'notspecified' or not license_url) else license_url
        return [LegalResource(id=url, title=title)]

    def _get_landing_page(self, dataset_dict):
        """Extracts the landing page URL if valid."""
        url = dataset_dict.get('url')
        return [Document(id=url)] if url and "https://" in str(url) else []

    def _get_measurement_technique(self, dataset_dict):
        """Returns a tuple (iri, label) for the measurement technique."""
        raw_iri = dataset_dict.get("measurement_technique_iri")
        raw_label = dataset_dict.get("measurement_technique")
        return (raw_iri or "http://purl.obolibrary.org/obo/OBI_0000070", raw_label or "assay")

    def _get_dates(self, dataset_dict):
        """Extracts and formats dates (YYYY-MM-DD)."""

        def clean_date(d):
            if not d:
                return None
            try:
                return str(d).split('T')[0]
            except Exception:
                return None

        return clean_date(dataset_dict.get('metadata_created')), clean_date(dataset_dict.get('metadata_modified'))

    def _get_theme(self, dataset_dict):
        """Returns the SKOS Concept LinkML object for the dcat:theme property, which for NFDI4Chem is always an instance
        from the EU Dataset Theme Vocabulary (http://publications.europa.eu/resource/authority/data-theme/) that is
        mandated by the DCAT-AP specification.
        TODO: We'll have to update DCAT-AP+ to allow providing the IRI of the concept directly instead of it as title.
        """
        theme = Concept(preferred_label="Science and technology",
                        title="http://publications.europa.eu/resource/authority/data-theme/TECH")
        return theme

class DCATAPPlusProfile(Helpers, EuropeanDCATAPProfile):
    """
    An RDF profile extending DCAT-AP for NFDI4Chem that inherits helper methods from NFDI4ChemHelpers.

    Extends the EuropeanDCATAPProfile to support NFDI4Chem-specific fields.
    """

    def parse_dataset(self, dataset_dict, dataset_ref):
        # TODO: Create a parser
        log.debug('Parsing dataset for NFDI4Chem')
        try:
            dataset_dict['title'] = str(dataset_ref.value(DCT.title))
            dataset_dict['notes'] = str(dataset_ref.value(DCT.description))
            dataset_dict['doi'] = str(dataset_ref.value(DCT.identifier))
            dataset_dict['language'] = [
                str(theme.value(SKOS.prefLabel)) for theme in dataset_ref.objects(DCAT.theme)
            ]
        except Exception as e:
            log.error(f"Error parsing dataset: {e}")
        return dataset_dict

    def graph_from_dataset(self, dataset_dict, dataset_ref):
        """
        Generates the RDF Graph for a dataset using DCAT-AP+ classes.
        """

        # 1. Bind Prefixes
        # Question from Philip to Bhavin: why do we need this here?
        # So far we only use the prefix map passed to the RDFLibDumper
        for prefix, namespace in namespaces.items():
            self.g.bind(prefix, namespace)

        # 2. Get Core IDs using Helpers
        dataset_id = self._get_dataset_id(dataset_dict)
        compound_id = self._get_compound_id(dataset_dict, dataset_id)
        sample_id = f"{dataset_id}#sample"
        meas_id = f"{dataset_id}#measurement"

        # 3. Load Schema (Cached, using Helper with DCAT-specific args)
        sv = self._get_schema_view(
            schema_name="dcat_ap_plus",
            local_filename="dcat_ap_plus.yaml",
            purl="https://w3id.org/nfdi-de/dcat-ap-plus/"
        )

        if not sv:
            log.critical("Cannot generate RDF: DCAT-AP+ Schema could not be loaded.")
            return

        # 4. Build Compound Entity
        compound = Entity(id=compound_id)

        # Qualitative Attributes
        if dataset_dict.get("inchi_key"):
            compound.has_qualitative_attribute.append(QualitativeAttribute(
                rdf_type=DefinedTerm(id='http://semanticscience.org/resource/CHEMINF_000059', title='InChiKey'),
                title="assigned InChIKey",
                value=dataset_dict.get("inchi_key")
            ))
        if dataset_dict.get("inchi"):
            compound.has_qualitative_attribute.append(QualitativeAttribute(
                rdf_type=DefinedTerm(id='http://semanticscience.org/resource/CHEMINF_000113', title='InChi'),
                title="assigned InChI",
                value=dataset_dict.get("inchi")
            ))
        if dataset_dict.get("smiles"):
            compound.has_qualitative_attribute.append(QualitativeAttribute(
                rdf_type=DefinedTerm(id='http://semanticscience.org/resource/CHEMINF_000018', title='SMILES'),
                title="assigned SMILES",
                value=dataset_dict.get("smiles")
            ))
        if dataset_dict.get("mol_formula"):
            compound.has_qualitative_attribute.append(QualitativeAttribute(
                rdf_type=DefinedTerm(id='http://semanticscience.org/resource/CHEMINF_000042',
                                     title='molecular formula'),
                title="assigned molecular formula",
                value=dataset_dict.get("mol_formula")
            ))
        if dataset_dict.get("iupacName"):
            compound.has_qualitative_attribute.append(QualitativeAttribute(
                rdf_type=DefinedTerm(id='http://semanticscience.org/resource/CHEMINF_000107', title='IUPAC name'),
                title="assigned IUPAC name",
                value=dataset_dict.get("iupacName")
            ))

        # Quantitative Attribute (Exact Mass) - Cast to Float
        if dataset_dict.get("exactmass"):
            try:
                mass_val = float(dataset_dict.get("exactmass"))
                compound.has_quantitative_attribute.append(QuantitativeAttribute(
                    rdf_type=DefinedTerm(id='http://semanticscience.org/resource/CHEMINF_000217',
                                         title='exact mass descriptor'),
                    has_quantity_type="http://qudt.org/vocab/quantitykind/MolarMass",
                    unit="https://qudt.org/vocab/unit/GM-PER-MOL",
                    title="exact mass",
                    value=mass_val
                ))
            except (ValueError, TypeError):
                log.warning(f"Invalid exactmass value '{dataset_dict.get('exactmass')}', skipping.")

        # 5. Build Sample Entity
        sample = EvaluatedEntity(
            id=sample_id,
            title='evaluated sample',
            has_part=[compound.id]
        )

        # 6. Build Measurement Activity
        tech_iri, tech_label = self._get_measurement_technique(dataset_dict)
        measurement = DataGeneratingActivity(
            id=meas_id,
            rdf_type=DefinedTerm(id=tech_iri, title=tech_label),
            evaluated_entity=[sample.id]
        )

        # 7. Build Dataset Object using Helpers
        creators = self._get_authors(dataset_dict)
        publisher = self._get_publisher(dataset_dict)
        legislation = self._get_license(dataset_dict, dataset_id)
        language = self._get_language(dataset_dict)
        landing_pages = self._get_landing_page(dataset_dict)
        release_date, mod_date = self._get_dates(dataset_dict)

        dataset = Dataset(
            id=dataset_id,
            title=dataset_dict.get("title"),
            description=self._get_description(dataset_dict),
            theme=self._get_theme(dataset_dict),
            identifier=dataset_id,
            other_identifier=self._get_other_ids(dataset_dict),
            release_date=release_date,
            modification_date=mod_date,
            creator=creators,
            language=[language],
            publisher=publisher,
            applicable_legislation=legislation,
            landing_page=landing_pages,
            conforms_to=Standard(title='DCAT-AP-PLUS', description='https://w3id.org/nfdi-de/dcat-ap-plus'),
            was_generated_by=[measurement.id],
            is_about_entity=[sample.id],
        )

        # 8. Serialize to Graph
        rdf_dumper = RDFLibDumper()
        prefix_map = {
            'CHEMINF': 'http://semanticscience.org/resource/CHEMINF_',
            'CHMO': 'http://purl.obolibrary.org/obo/CHMO_',
            'CHEBI': 'http://purl.obolibrary.org/obo/CHEBI_'
        }

        try:
            graph = rdf_dumper.as_rdf_graph(dataset, schemaview=sv, prefix_map=prefix_map)
            graph += rdf_dumper.as_rdf_graph(sample, schemaview=sv, prefix_map=prefix_map)
            graph += rdf_dumper.as_rdf_graph(compound, schemaview=sv, prefix_map=prefix_map)
            graph += rdf_dumper.as_rdf_graph(measurement, schemaview=sv, prefix_map=prefix_map)

            for triple in graph:
                self.g.add(triple)
        except Exception as e:
            log.error(f"RDF Serialization failed: {e}")

