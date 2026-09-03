# -*- coding: utf-8 -*-
"""Compatibility helpers for datetime values rendered by CKAN 2.9."""

import datetime
import re


_ISO_DATETIME = re.compile(
    r'^\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2}'
    r'(?:\.\d{1,6})?(?:Z|[+-]\d{2}:\d{2})?$'
)


def datetime_for_ckan_form(value):
    """Return a naive UTC ISO value accepted by CKAN 2.9, or ``None``.

    Scheming accepts ISO-8601 timezone suffixes when saving, while CKAN 2.9's
    ``date_str_to_datetime`` (used by its edit form) does not.  Normalize only
    for rendering; the stored metadata is not modified.
    """
    if value is None or value == '':
        return value

    if isinstance(value, datetime.datetime):
        parsed = value
    else:
        if not isinstance(value, str) or not _ISO_DATETIME.match(value):
            return None
        parse_value = value[:-1] + '+00:00' if value.endswith('Z') else value
        has_tz = value.endswith('Z') or bool(re.search(r'[+-]\d{2}:\d{2}$', value))
        try:
            parsed = datetime.datetime.strptime(
                parse_value,
                '%Y-%m-%dT%H:%M:%S' +
                ('.%f' if '.' in parse_value else '') +
                ('%z' if has_tz else ''),
            )
        except (TypeError, ValueError, OverflowError):
            try:
                parsed = datetime.datetime.strptime(
                    parse_value,
                    '%Y-%m-%d %H:%M:%S.%f' if '.' in parse_value
                    else '%Y-%m-%d %H:%M:%S',
                )
            except (TypeError, ValueError, OverflowError):
                return None

    if parsed.tzinfo is not None:
        parsed = parsed.astimezone(datetime.timezone.utc).replace(tzinfo=None)

    return parsed.isoformat()
