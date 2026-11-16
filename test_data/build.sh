#!/usr/bin/env bash
# Rebuilds the walk-file fixtures from their GFAs with gfa-to-tabix 5.0, in
# walk format 8, and the fixtures of extension/ with that repo's corrupt.py.
# Usage: test_data/build.sh <gfa-to-tabix binary> <gfa-to-tabix checkout>
set -euo pipefail
cd "$(dirname "$0")"
B=${1:?the gfa-to-tabix binary}
CORRUPT=${2:?the gfa-to-tabix checkout}/docs/departures/scripts/broken/corrupt.py
build() {
  local name=$1 gfa=$2; shift 2
  "$B" "$gfa" --walks -o "$name.tmp" "$@" > /dev/null 2>&1
  local out
  for out in "$name".tmp.*.walks.gz; do
    local sample=${out#"$name".tmp.}; sample=${sample%.walks.gz}
    mv "$out" "$name.$sample.walks.gz"
    mv "$out.idx" "$name.$sample.walks.gz.idx"
  done
}
build header header.GRCh38.gfa --refs GRCh38 --tile 1000
build maxnode maxnode.GRCh38.gfa --refs GRCh38 --tile 1000
build revisit revisit.GRCh38.gfa --refs GRCh38 --tile 1000
build farpass farpass.GRCh38.gfa --refs GRCh38 --tile 1000
build follow follow.GRCh38.gfa --refs GRCh38 --tile 1000
build rejoin follow.GRCh38.gfa --refs GRCh38 --tile 10000
build records records.GRCh38.gfa --refs GRCh38 --tile 1000
build skip skip.GRCh38.gfa --refs GRCh38 --tile 10
build span span.gfa --refs GRCh38 --tile 20
build attach attach.GRCh38.gfa --refs GRCh38 --tile 20
build unabridged unabridged.R.gfa --refs R --tile 10
build abridged unabridged.R.gfa --refs R --tile 10 --abridge-bytes 1
build interleaved interleaved.R.gfa --refs R --tile 10 --stay 0 --abridge-steps 3
build ownends ownends.R.gfa --refs R --tile 10 --abridge-bytes 1
build contigends contigends.R.gfa --refs R --tile 10 --abridge-bytes 1
build hprc-chrM hprc-chrM.GRCh38.gfa --refs GRCh38 --sequences
build jumpback jumpback.GRCh38.gfa --refs GRCh38 --tile 20
build loop loop.R.gfa --refs R --tile 10
build ownzero ownzero.R.gfa --refs R --tile 10
build excursion excursion.R.gfa --refs R --tile 64
build hole hole.R.gfa --refs R --tile 16
build holetile hole.R.gfa --refs R --tile 64
build refhole refhole.R.gfa --refs R --tile 64
build abut abut.R.gfa --refs R --tile 20
build copies copies.R.gfa --refs R --tile 1000
build copyends copyends.R.gfa --refs R --tile 1000 --abridge-bytes 1
build zerostep zerostep.R.gfa --refs R --tile 8
build early early.R.gfa --refs R --tile 1000
build twolines twolines.R.gfa --refs R --tile 20
build abredge abredge.R.gfa --refs R --tile 20 --abridge-bytes 1
build drawnedge drawnedge.R.gfa --refs R --tile 10000 --abridge-steps 3
build attachmany attachmany.R.gfa --refs R --tile 10
build blocks blocks.R.gfa --refs R --tile 60 --abridge-bytes 1
build between between.R.gfa --refs R --tile 60 --abridge-bytes 1
build unbuilt unbuilt.R.gfa --refs R --tile 60 --abridge-bytes 1
build hairpin hairpin.R.gfa --refs R --tile 36 --abridge-bytes 1
build revref revref.gfa --refs GRCh38 --abridge-bytes 1
build cover cover.R.gfa --refs R --abridge-bytes 1
build counts counts.R.gfa --refs R --tile 20
build countsabridged counts.R.gfa --refs R --tile 20 --abridge-steps 4
build twocontigs twocontigs.R.gfa --refs R --tile 20
build contained contained.R.gfa --refs R --tile 40 --abridge-bytes 1
for twin in twinmidlineA twinmidlineB twinotherseqA twinotherseqB twinown690 twinown700 pickwrong twopieces twopiecesonetile; do
  build "$twin" "$twin.R.gfa" --refs R --tile 12 --stay 0 --abridge-bytes 1
done
build adv adv.R.gfa --refs R --tile 60 --abridge-bytes 1
build edge edge.R.gfa --refs R --tile 60 --abridge-bytes 1
build farrejoin farrejoin.R.gfa --refs R --tile 10
build detailtwice detailtwice.R.gfa --refs R --tile 40 --abridge-bytes 1
# The older twins, unabridged and the *whole files without --abridge-steps 0
# --abridge-bytes 0, keep the default abridge budgets, which abridge no tile of
# a graph this small, so each holds every walk as one built with both 0 does
build detailtwicewhole detailtwice.R.gfa --refs R --tile 40
build sortsfirst sortsfirst.R.gfa --refs R --tile 40 --abridge-bytes 1
build sortsfirstwhole sortsfirst.R.gfa --refs R --tile 40
build ownendswhole ownends.R.gfa --refs R --tile 10
build unbuiltwhole unbuilt.R.gfa --refs R --tile 60
build edgewhole edge.R.gfa --refs R --tile 60
build advwhole adv.R.gfa --refs R --tile 60
build gapnan gapnan.R.gfa --refs R --tile 20 --abridge-bytes 1
build gapnanwhole gapnan.R.gfa --refs R --tile 20 --abridge-steps 0 --abridge-bytes 0
build secondline secondline.R.gfa --refs R --tile 20 --abridge-bytes 1
build secondlinewhole secondline.R.gfa --refs R --tile 20 --abridge-steps 0 --abridge-bytes 0
build extension/plain extension/follow.R.gfa --refs R --tile 100 --sequences --abridge-bytes 1
addition() {
  local name=$1; shift
  python3 -I "$CORRUPT" extension/plain.R.walks.gz "extension/$name.R.walks.gz" "$@" > /dev/null
}
addition tags --tag zz:i:1
addition unit-kind $'#unit\tR#0#c\t0\t36\ttable\t.\n' $'#unit\tR#0#c\t0\t36\tfuture\tx\nR#0#c\t0\t36\tfuture\tanything\n#unit\tR#0#c\t0\t36\ttable\t.\n'
addition header-line $'#reference\tR\n' $'#reference\tR\n#future\tx\n'
