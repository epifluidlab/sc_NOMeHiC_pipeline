my $p=$ARGV[0] or die;
my $input=$ARGV[1] or die;
my $genome=$ARGV[2] || "/projects/b1198/epifluidlab/eric/projects/scnomehic_opt/reference/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.chrom.sizes";
my $methy="$p.methy.bedgraph";
my $cov="$p.cov.bedgraph";
my $m_count="$p.methy_count.bedgraph";


	if($input=~/\.gz$/){
		`zcat $input | awk '\$7!="."'| grep -v "_"| cut -f1-3,7-8 | grep -v \"^track\" | perl -ne 'chomp;\@f=split \"\\t\";\$count=int(\$f[3]*\$f[4]/100);print \"\$f[0]\\t\$f[1]\\t\$f[2]\\t\$count\\n\";' > $m_count`;
		`zcat $input | awk '\$7!="."'| grep -v "_"| cut -f1-3,8 | grep -v \"^track\" > $cov`;
		`zcat $input | awk '\$7!="."'| grep -v "_"| cut -f1-3,7 | grep -v \"^track\" > $methy`;
	}else{
		`cut -f1-3,7-8 $input | grep -v "_"| awk '\$4!="."' | grep -v \"^track\" | perl -ne 'chomp;\@f=split \"\\t\";\$count=int(\$f[3]*\$f[4]/100);print \"\$f[0]\\t\$f[1]\\t\$f[2]\\t\$count\\n\";' > $m_count`;
		`cut -f1-3,8 $input | grep -v "_"| awk '\$4!=0'| grep -v \"^track\" > $cov`;
		`cut -f1-3,7 $input | grep -v "_" | awk '\$4!="."'| grep -v \"^track\" > $methy`;

	}

my $methy_bw="$p.methy.bw";
my $cov_bw="$p.cov.bw";
my $m_count_bw="$p.methy_count.bw";

`bedGraphToBigWig $methy $genome $methy_bw`;
`bedGraphToBigWig $cov $genome $cov_bw`;
`bedGraphToBigWig $m_count $genome $m_count_bw`;

`unlink $methy`;
`unlink $cov`;
`unlink $m_count`;



