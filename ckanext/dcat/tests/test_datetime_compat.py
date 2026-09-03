# -*- coding: utf-8 -*-

import pytest

from ckanext.dcat.datetime_compat import datetime_for_ckan_form


@pytest.mark.parametrize('value, expected', [
    ('2026-09-03T10:09:43Z', '2026-09-03T10:09:43'),
    ('2026-09-03T10:09:43.123Z', '2026-09-03T10:09:43.123000'),
    ('2026-09-03T10:09:43+00:00', '2026-09-03T10:09:43'),
    ('2026-09-03T10:09:43.123+00:00', '2026-09-03T10:09:43.123000'),
    ('2026-09-03T12:09:43+02:00', '2026-09-03T10:09:43'),
    ('2026-09-03 10:09:43', '2026-09-03T10:09:43'),
    ('', ''),
    (None, None),
    ('not-a-date', None),
])
def test_datetime_for_ckan_form(value, expected):
    assert datetime_for_ckan_form(value) == expected
