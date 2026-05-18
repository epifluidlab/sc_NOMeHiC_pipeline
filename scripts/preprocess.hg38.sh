zcat -f -- "$1" | cut -d' ' -f1-7 | sed 's/ /\t/g'| perl -ne 'chomp;@f=split "\t";print join("\t",@f)."\n";' | gzip -c > "$2"
