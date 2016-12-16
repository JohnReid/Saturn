#!/usr/bin/env python2

import pprint
from BCBio.GFF import GFFExaminer, parse
import gzip

in_file = "TSS_human_with_gencodetss_notlow_ext50eachside_merged_withgenctsscoord_andgnlist.gff.gz"
examiner = GFFExaminer()
with gzip.open(in_file) as in_handle:
    for rec in parse(in_handle):
        print rec
    # pprint.pprint(examiner.available_limits(in_handle))
    # pprint.pprint(examiner.parent_child_map(in_handle))
