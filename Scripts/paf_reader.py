#!usr/bin/env python3

"""
Script for reworking rhe Paf file.
"""

from __future__ import division
from collections import namedtuple


__all__ = ["parse_paf"]

try:
    import pandas as pd
except Exception as E:
    pandas = False
    e = E
else:
    pandas = True


class _PAF:
    """Base PAF methods, can't guarantee field names here so use indices"""

    def __str__(self):
        """Formats a record as a PAF line for writing to a file"""
        return "{}\t{}".format("\t".join(map(str, self[:-1])), self._fmt_tags())

    def _fmt_tags(self):
        """Format tag dict as SAM style"""
        return "\t".join("{}:{}:{}".format(*t) for t in self[-1].values())

    def blast_identity(self):
        """BLAST identity, see:
        https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
        """
        return self[9] / self[10]


SAM_TAG = namedtuple("tag", ["name", "type", "value"])
FIELDS = ["query_name", "query_length", "query_start", "query_end",
          "strand", "target_name", "target_length", "target_start",
          "target_end", "residue_matches", "alignment_block_length",
          "mapping_quality", "tags"]
NA_VALUES = ["*"]
SAM_TYPES = {"i": int, "A": str, "f": float, "Z": str}


def _expand_dict_in_series(df, field):
    """Convert a Series of dict to Series and add to the original DataFrame
    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame with a Series of dict
    field : str
        The Series of dicts to expand
    Returns
    -------
    pd.DataFrame
        The orignal DataFrame with extra Series from the dicts
    """
    return df.join(
        pd.DataFrame(
            [{k: v for k, _, v in r.values()} for r in df.pop(field).tolist()]
        ),
        rsuffix="_tag",
    )


def _parse_tags(tags):
    """ Convert a list of SAM style tags, from a PAF file, to a dict """
    return {
        tag: SAM_TAG(tag, type_, SAM_TYPES.get(type_, lambda x: x)(val))
        for tag, type_, val in (x.split(":", 2) for x in tags)
    }


def _paf_generator(file_like, fields=None, na_values=None, na_rep=None):
    """ Generator that returns namedtuples from a PAF file """
    if len(fields) != 13:
        raise ValueError("{} fields provided, expected 13".format(len(fields)))
    _PAF_nt = namedtuple("PAF", fields)
    PAF = type("PAF", (_PAF, _PAF_nt), dict())
    for record in file_like:
        record = record.strip()
        if not record:
            continue
        record = record.split("\t")
        yield PAF(
            str(record[0]),
            int(record[1]) if record[1] not in na_values else na_rep,
            int(record[2]) if record[2] not in na_values else na_rep,
            int(record[3]) if record[3] not in na_values else na_rep,
            str(record[4]),
            str(record[5]),
            int(record[6]) if record[6] not in na_values else na_rep,
            int(record[7]) if record[7] not in na_values else na_rep,
            int(record[8]) if record[8] not in na_values else na_rep,
            int(record[9]) if record[9] not in na_values else na_rep,
            int(record[10]) if record[10] not in na_values else na_rep,
            int(record[11]) if record[11] not in na_values else na_rep,
            _parse_tags(record[12:]),
        )


def parse_paf(file_like, fields=None, na_values=None, na_rep=0, dataframe=False):
    """ Read a minimap2 PAF file as either an iterator or a pandas.DataFrame """
    fields = FIELDS if fields is None else fields
    na_values = set(NA_VALUES if na_values is None else na_values + NA_VALUES)
    if not isinstance(na_rep, (int, float)):
        raise ValueError("na_rep must be int or float")

    if dataframe and pandas:
        df = pd.DataFrame(
            (line.strip().split("\t", 12) for line in file_like if line.strip()),
            columns=fields,
        )
        df = df.join(
            pd.DataFrame(
                df.pop(fields[-1])
                .str.findall(r"([^\t]+?):[A-Za-z]+?:(.+?)")
                .map(dict)
                .to_list()
            ),
            rsuffix="_tag",
        )
        if df.empty:
            return pd.DataFrame(columns=fields)
        df = df.replace(
            {
                fields[i]: {v: na_rep for v in na_values}
                for i in (2, 3, 4, 7, 8, 9, 10, 11, 12)
            }
        )
        return df.infer_objects()
    elif dataframe and not pandas:
        raise ImportError(e)
    else:
        return _paf_generator(file_like, fields, na_values, na_rep)
