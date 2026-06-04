import logging
from rdflib import term, URIRef, BNode, Literal, Graph
import ckantoolkit as toolkit
from .base import (RDF,
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
                   CHEMINF,
                   CHMO,
                   OBI,
                   IAO,
                   PROV,
                   CHEBI,
                   NMR,
                   QUDT,
                   NCIT,
                   FIX,
                   namespaces,
                   )
from linkml_runtime.dumpers import RDFLibDumper

# Import base class - the DCAT profile we are inheriting from
from .dcat_ap_plus import DCATAPPlusProfile
from .euro_dcat_ap import EuropeanDCATAPProfile

# Import ChemDCAT-AP specific dataclasses (Local copy for Python 3.7 compatibility)
# NOTE: In the future, replace this with: from chem_dcat_ap.datamodel.chem_dcat_ap import ...
from .dcat_4c_ap import (
    SubstanceSample,
    SubstanceSampleCharacterizationDataset,
    SubstanceSampleCharacterization,
    InChi, InChIKey, IUPACName, SMILES, MolecularFormula, MolarMass,
    ChemicalEntity, DefinedTerm, Standard
)

log = logging.getLogger(__name__)


class ChemDCATAPProfile(DCATAPPlusProfile):
    """
    ChemDCAT-AP Profile.
    Inherits all data extraction and helper logic from DCATAPPlusProfile.
    Only implements the specific ChemDCAT-AP graph construction.
    """
    def parse_dataset(self, dataset_dict, dataset_ref):
        """
        Parses the RDF reference back into a dictionary.
        Re-using logic compatible with the parent class.
        TODO: Create a parser
        """
        log.debug("Parsing dataset for ChemDCAT-AP")
        try:
            dataset_dict["title"] = str(dataset_ref.value(DCT.title) or "")
            dataset_dict["notes"] = str(dataset_ref.value(DCT.description) or "")
            dataset_dict["doi"] = str(dataset_ref.value(DCT.identifier) or "")
            # Handle language list or single value gracefully
            lang_objs = list(dataset_ref.objects(DCT.language))
            if lang_objs:
                dataset_dict["language"] = str(lang_objs[0].value(SKOS.prefLabel) or lang_objs[0])
            else:
                dataset_dict["language"] = ""
        except Exception as e:
            log.error(f"Error parsing dataset: {e}")
        return dataset_dict


    def graph_from_dataset(self, dataset_dict, dataset_ref):
        """
        Generates the RDF Graph for a dataset using ChemDCAT-AP classes.
        """

        # 1. Bind Prefixes
        # Question from Philip to Bhavin: why do we need this here?
        # So far we only use the prefix map passed to the RDFLibDumper
        for prefix, namespace in namespaces.items():
            self.g.bind(prefix, namespace)

        # 2. Get Core IDs using Inherited Helpers
        dataset_id = self._get_dataset_id(dataset_dict)
        compound_id = self._get_compound_id(dataset_dict, dataset_id)
        sample_id = f"{dataset_id}#sample"
        meas_id = f"{dataset_id}#measurement"

        # 3. Load Schema (Cached, using Inherited Helper with Chem-specific args)
        sv = self._get_schema_view(
            schema_name="chem_dcat_ap",
            local_filename="chem_dcat_ap.yaml", # Fallback for now, until PURL is live
            purl="https://w3id.org/nfdi-de/dcat-ap-plus/chemistry/"
        )

        if not sv:
            log.critical("Cannot generate RDF: ChemDCAT-AP Schema could not be loaded.")
            return

        # 4. Build ChemicalEntity (Chem-Specific)
        compound_kwargs = {"id": compound_id}

        if dataset_dict.get("inchi_key"):
            compound_kwargs["inchikey"] = InChIKey(
                title="assigned InChIKey",
                value=dataset_dict.get("inchi_key")
            )
        if dataset_dict.get("inchi"):
            compound_kwargs["inchi"] = InChi(
                title="assigned InChI",
                value=dataset_dict.get("inchi")
            )
        if dataset_dict.get("smiles"):
            compound_kwargs["smiles"] = SMILES(
                title="assigned SMILES",
                value=dataset_dict.get("smiles")
            )
        if dataset_dict.get("mol_formula"):
            compound_kwargs["molecular_formula"] = MolecularFormula(
                title="assigned molecular formula",
                value=dataset_dict.get("mol_formula")
            )
        if dataset_dict.get("exactmass"):
            # Cast to float for type correctness (matches script & DCAT profile)
            try:
                mass_val = float(dataset_dict.get("exactmass"))
                compound_kwargs["has_molar_mass"] = MolarMass(
                    has_quantity_type="http://qudt.org/vocab/quantitykind/MolarMass",
                    unit="https://qudt.org/vocab/unit/GM-PER-MOL",
                    title="assigned exact mass",
                    value=mass_val
                )
            except (ValueError, TypeError):
                log.warning(f"Invalid exactmass value '{dataset_dict.get('exactmass')}', skipping.")

        if dataset_dict.get("iupacName"):
            compound_kwargs["iupac_name"] = IUPACName(
                title="assigned IUPAC name",
                value=dataset_dict.get("iupacName")
            )

        compound_chem = ChemicalEntity(**compound_kwargs)

        # 5. Build SubstanceSample (Chem-Specific)
        sample_chem = SubstanceSample(
            id=sample_id,
            title="evaluated sample",
            composed_of=[compound_chem.id]
        )

        # 6. Build Measurement (Chem-Specific)
        tech_iri, tech_label = self._get_measurement_technique(dataset_dict)
        measurement_chem = SubstanceSampleCharacterization(
            id=meas_id,
            description="The kind of activity/process used to generate the dataset",
            rdf_type=DefinedTerm(id=tech_iri, title=tech_label),
            evaluated_entity=[sample_chem.id]
        )

        # 7. Build Dataset (Chem-Specific Class, Shared Metadata Helpers)
        creators = self._get_authors(dataset_dict)
        publisher = self._get_publisher(dataset_dict)
        legislation = self._get_license(dataset_dict, dataset_id)
        language = self._get_language(dataset_dict)
        landing_pages = self._get_landing_page(dataset_dict)
        release_date, mod_date = self._get_dates(dataset_dict)

        dataset_chem = SubstanceSampleCharacterizationDataset(
            id=dataset_id,
            title=dataset_dict.get("title"),
            description=self._get_description(dataset_dict),
            theme=self._get_theme(dataset_dict),
            identifier=dataset_id,
            other_identifier= self._get_other_ids(dataset_dict),
            release_date=release_date,
            modification_date=mod_date,
            creator=creators,
            language=[language],
            publisher=publisher,
            applicable_legislation=legislation,
            landing_page=landing_pages,
            conforms_to=Standard(title='ChemDCAT-AP',
                                 description='https://w3id.org/nfdi-de/dcat-ap-plus/chemistry/'),
            was_generated_by=[measurement_chem.id],
            is_about_entity=[sample_chem.id],
        )

        # 8. Serialize to Graph
        rdf_dumper = RDFLibDumper()
        prefix_map = {
            'CHEMINF': 'http://semanticscience.org/resource/CHEMINF_',
            'CHMO': 'http://purl.obolibrary.org/obo/CHMO_',
            'CHEBI': 'http://purl.obolibrary.org/obo/CHEBI_'
        }

        try:
            graph = rdf_dumper.as_rdf_graph(dataset_chem, schemaview=sv, prefix_map=prefix_map)
            graph += rdf_dumper.as_rdf_graph(sample_chem, schemaview=sv, prefix_map=prefix_map)
            graph += rdf_dumper.as_rdf_graph(compound_chem, schemaview=sv, prefix_map=prefix_map)
            graph += rdf_dumper.as_rdf_graph(measurement_chem, schemaview=sv, prefix_map=prefix_map)

            for triple in graph:
                self.g.add(triple)
        except Exception as e:
            log.error(f"ChemDCAT-AP RDF Serialization failed: {e}")