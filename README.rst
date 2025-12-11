========
FeaVar
========

.. image:: https://img.shields.io/pypi/v/FeaVar.svg
        :target: https://pypi.python.org/pypi/FeaVar

.. image:: https://readthedocs.org/projects/FeaVar/badge/?version=latest
        :target: https://FeaVar.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

FeaVar is a python package that computes clusters of unique substrings or "variant types" (SFVTs) for user-selected sequence features in a set of aligned sequences.

This package is a generalization of the sequence feature variant type (SFVT) work undertaken in the laboratory of Dr. Richard Scheuermann while at UT Southwestern.

Features
--------

* **SFVT Computation:** Identifies unique Sequence Feature Variant Types (SFVTs) based on user-defined substrings (e.g., specific genetic motifs or domains).
* **Reference Correction:** Automatically adjusts position coordinates to account for gaps and insertions relative to a reference sequence in the alignment.
* **Metadata Integration:** Merges computed variant types with external metadata (CSV/TSV) to analyze relationships between variants and sample attributes.
* **Visualization:** Automatically generates bar charts and stacked plots for the top occurring variants.
* **Flexible Input:** Supports Clustal and other alignment formats via Biopython.
* **Exportable Results:** Outputs all processed data, including variant counts and merged metadata, to CSV files for downstream analysis.

Installation
------------

To install FeaVar locally, run the following command from the root of the repository:

.. code-block:: bash

    pip install .

Usage
-----

Once installed, you can run ``feavar`` from the command line.

**Basic Syntax**

.. code-block:: bash

    feavar -a <ALIGNMENT_FILE> -r <REF_ID> -p <POSITIONS>

**Arguments**

* ``-a``, ``--alignment``: Path to the multiple sequence alignment file (required).
* ``-r``, ``--reference_identifier``: The ID of the reference sequence as it appears in the alignment (required).
* ``-p``, ``--positions``: The positions to analyze (1-based). Can be a range (``100-110``) or comma-separated list (``100-110,120,130``).
* ``-m``, ``--metadata_file``: Path to a tab-delimited metadata file to merge with results (optional).
* ``-t``, ``--top``: Number of top variant types to plot (default: 10).
* ``-f``, ``--alignment_format``: Format of the alignment file (default: ``clustal``).

**Example**

.. code-block:: bash

    feavar -a data/flu_sequences.clw -r CY021716 -p "124-142" -m data/metadata.tsv -t 15

References
----------

- Noronha, J. M., Liu, M., Squires, R. B., Pickett, B. E., Hale, B. G., Air, G. M., et al. (2012). Influenza virus sequence feature variant type analysis: evidence of a role for NS1 in influenza virus host range restriction. *Journal of Virology*, 86(10), 5857–5866. http://doi.org/10.1128/JVI.06901-11
- Karp, D. R., Marthandan, N., Marsh, S. G. E., Ahn, C., Arnett, F. C., DeLuca, D. S., et al. (2009). Novel sequence feature variant type analysis of the HLA genetic association in systemic sclerosis. *Human Molecular Genetics*, 19(4), 707–719. http://doi.org/10.1093/hmg/ddp521
- Thomson, Glenys & Marthandan, Nishanth & Hollenbach, Jill & Mack, St & A Erlich, Henry & Single, Richard & Waller, Matthew & Marsh, Steven & A Guidry, Paula & R Karp, David & Scheuermann, Richard & Thompson, Susan & N Glass, David & Helmberg, Wolfgang. (2010). Sequence feature variant type (SFVT) analysis of the HLA genetic association in juvenile idiopathic arthritis. *Pacific Symposium on Biocomputing*. 359-70. 10.1142/9789814295291_0038

License
-------

This project is licensed under the MIT License.