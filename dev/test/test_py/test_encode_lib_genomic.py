"""Make sure that
    - pipeline's src/ directory is added to PYTHONPATH.
    - pipeline's conda environment is activated.
        Or run this in pipeline's docker container.

This module does not cover all functions defined in encode_lib_genomic.
It only covers the following newly added functions.
"""

import gzip
import pytest
from encode_lib_genomic import bed_clip
from textwrap import dedent

CHRSZ_HG38 = dedent("""\
    chr1\t248956422
    chr2\t242193529
    chr3\t198295559
    chr4\t190214555
    chr5\t181538259
    chr6\t170805979
    chr7\t159345973
    chr8\t145138636
    chr9\t138394717
    chr10\t133797422
    chr11\t135086622
    chr12\t133275309
    chr13\t114364328
    chr14\t107043718
    chr15\t101991189
    chr16\t90338345
    chr17\t83257441
    chr18\t80373285
    chr19\t58617616
    chr20\t64444167
    chr21\t46709983
    chr22\t50818468
    chrX\t156040895
    chrY\t57227415
    chrM\t16569
""")


def test_bed_clip(tmp_path):
    chrsz = tmp_path / 'chrsz'
    chrsz.write_text(CHRSZ_HG38)

    out_bed = tmp_path / 'out.bed'

    # no change if it's inside 0-chrsz
    bed_in = tmp_path / 'test.in.bed'
    bed_in.write_text('chrM\t2\t10000\n')
    bed_clip(str(bed_in), str(chrsz), str(out_bed), no_gz=True)
    assert out_bed.read_text() == 'chrM\t2\t10000\n'

    # out-of-bound error should occur
    bed_oob1 = tmp_path / 'test.oob1.bed'
    bed_oob1.write_text('chrM\t-10\t-1\n')
    with pytest.raises(Exception):
        bed_clip(str(bed_oob1), str(chrsz), str(out_bed), no_gz=True)

    # should be truncated to 0 (left) or chrsz (right)
    bed_left_crossing = tmp_path / 'test.lc.bed'
    bed_left_crossing.write_text('chrM\t-1\t13000\n')
    bed_right_crossing = tmp_path / 'test.rc.bed'
    bed_right_crossing.write_text('chrM\t13000\t17000\n')
    bed_larger = tmp_path / 'test.l.bed'
    bed_larger.write_text('chrM\t-1\t17000\n')    

    bed_clip(str(bed_left_crossing), str(chrsz), str(out_bed), no_gz=True)
    assert out_bed.read_text() == 'chrM\t0\t13000\n'

    bed_clip(str(bed_right_crossing), str(chrsz), str(out_bed), no_gz=True)
    assert out_bed.read_text() == 'chrM\t13000\t16569\n'

    bed_clip(str(bed_larger), str(chrsz), str(out_bed), no_gz=True)
    assert out_bed.read_text() == 'chrM\t0\t16569\n'

    # test no_gz flag
    out_bed = tmp_path / 'out.bed.gz'
    bed_clip(str(bed_larger), str(chrsz), str(out_bed), no_gz=False)
    with gzip.open(str(out_bed), 'rb') as fp:
        assert fp.read().decode() == 'chrM\t0\t16569\n'

