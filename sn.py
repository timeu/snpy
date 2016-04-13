# -*- coding: utf-8 -*-
"""
    snpy
    ~~~~

    ``snpy`` provides a fancy API for accessing `openSNP`_ data. The
    current implementation only supports working with *local* or downloaded
    files, but JSON API interaction is planned.

    .. _openSNP: http://opensnp.org

    :copyright: 2012 by Sergei Lebedev
    :license: WTFPL, see LICENSE for details
"""

import csv
import os.path
from collections import namedtuple

try:
    import vcf
except ImportError:
    vcf = None  # 23andMe exome data won't be supported.


class SNPyError(Exception):
    """Generic ``snpy`` exception class."""


class UnknownSource(SNPyError):
    """Raised when a parsed file has unknown source type."""


_SNP = namedtuple("_SNP", ["name",
                           "variation",
                           "chromosome",
                           "position",
                           "strand",
                           "genotype"])


class SNP(_SNP):
    """A wrapper for SNP data, provided by various formats."""
    __slots__ = ()
    
    def __new__(cls, name, chromosome, position, genotype,
                variation=None, strand=None):
        return super(SNP, cls).__new__(cls, name, variation, chromosome,
                                       int(position), strand, genotype)


def _23andme(fd):
    handle = csv.DictReader(fd,
        fieldnames=["name", "chromosome", "position", "genotype"],
        delimiter="\t")

    for row in handle:
        if not row["name"].startswith("#"):
            yield SNP(**row)


def _23andme_exome(fd):
    if vcf is None:
        raise RuntimeError("PyVCF not available, please 'easy_install' it.")

    for r in vcf.VCFReader(fd):
        if not r.is_snp:
            continue  # XXX Is it even possible?

        for sample in r.samples:
            yield SNP(name=r.ID, chromosome=r.CHROM, position=r.POS,
                      genotype=sample.gt_bases.replace("/", ""))


def decodeme(fd):
    handle = csv.DictReader(fd,
        fieldnames=["name", "variation", "chromosome", "position",
                    "strand", "genotype"])

    for row in handle:
        # A flanky header criterion -- the last column should be
        # 'XX', where X is one of "ACGT-".
        if len(row["genotype"]) == 2:
            yield SNP(**row)


def ftdna(fd):
    handle = csv.DictReader(fd,
        fieldnames=["name", "chromosome", "position", "genotype"])

    for row in handle:
        if row["position"].isdigit():
            yield SNP(**row)

def ancestry(fd):
    handle = csv.DictReader(fd,
        fieldnames=["name", "chromosome", "position", "allele1","allele2"],
        delimiter="\t")
    for row in handle:
        if not row["name"].startswith("#") and row["position"].isdigit():
            row['genotype'] = row['allele1']+row['allele2']
            del(row['allele1'])
            del(row['allele2'])
            yield SNP(**row)
    
def guess_source(path):
    name, ext = os.path.split(path)
    if ext == "vcf":
        return ext  # VCF is easy ;)

    # Okay, maybe it's in openSNP format: ``format.submission_id``.
    try:
        source, _ = path.rsplit(os.path.extsep, 2)[-2:]
    except ValueError:
        raise UnknownSource(path)
    else:
        return source


def _guess_source_from_generator(it):
    first_line = next(it)
    # check if 23andMe or ancestry
    if first_line[0] == '#':
        if 'AncestryDNA' in first_line:
            return 'ancestry'
        elif '23andMe' in first_line:
            return '23andme'
        # if can't match first line iterate until header line and check there
        for row in it:
            if row.startswith("#"): 
                continue
            fields = row.split('\t')
            if len(fields) == 5:
                return 'ancestry'
            else:
                return '23andme'
    else:
        fields = first_line.split(',')
        if len(fields) == 6 and fields == ['Name','Variation','Chromosome','Position','Strand','YourCode']:
            return 'decodeme'
        elif len(fields) == 5 and fields == ['RSID','CHROMOSOME','POSITION','RESULT']:
            return 'ftdna'
    raise UnknownSource(path)
   
    

def guess_source_from_file(path):
    if not hasattr(path, 'read'):
        with open(path,'r') as fd:
            return  _guess_source_from_generator(fd)
    else: 
        return _guess_source_from_generator(path)

def guess_source_from_content(content):
    return _guess_source_from_generator(iter(content.splitlines()))


def parse(path, source=None):
    """Returns a generator yielding :class:`SNP` from an openSNP file
    at a given location.

    :param str path: path to openSNP file, *all* formats are supported.
    :param str source: should be one of ``"23andme"``, ``"vcf"``,
        ``"decodeme"`` or ``"ftdna"``; `source`` will be extracted from
        the filename, if not provided explicitly.
    :returns list: of :class:`SNP` instances from a given file.
    :raises RuntimeError: if a given file cannot be parsed.
    """
    if source is None:
        source = guess_source(path)
        
    try:
        handler = {"23andme": _23andme,
                   "23andme-exome-vcf": _23andme_exome,
                   "ftdna-illumina": ftdna,
                   "decodeme": decodeme,
                   "vcf": _23andme_exome,
                   "ftdna": ftdna,
                   "ancestry":ancestry}[source]
    except KeyError:
        raise UnknownSource(path)
    else:
        if not hasattr(path, 'read'):
            with open(path,'r') as fd:
                return handler(fd)
        else: 
            return handler(path)
     
        
