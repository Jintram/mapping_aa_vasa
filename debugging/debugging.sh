NS500414:716:HTMH2BGXF:1:11202:1204:11855;SS:TACTGGTA


# problem lines are (as they look in bed file):
14	23426826	23426863	NS500414:716:HTMH2BGXF:1:11202:1204:11855;SS:TACTGGTA;CB:TACTGGTA;QT:eeeVeeee;RX:CGAGAG;RQ:aaaaae;SM:011;CL:targeted;p_seq:AACGTACAAAGTGGGGATGG;C1:37M;C2:37M375N36M;nM:0/1	255	+
14	23426827	23427275	NS500414:716:HTMH2BGXF:1:11202:1204:11855;SS:TACTGGTA;CB:TACTGGTA;QT:eeeVeeee;RX:CGAGAG;RQ:aaaaae;SM:011;CL:targeted;p_seq:AACGTACAAAGTGGGGATGG;C1:37M;C2:37M375N36M;nM:0/2	255	-


# let's test my awk command on this:
myline="14	23426826	23426863	NS500414:716:HTMH2BGXF:1:11202:1204:11855;SS:TACTGGTA;CB:TACTGGTA;QT:eeeVeeee;RX:CGAGAG;RQ:aaaaae;SM:011;CL:targeted;p_seq:AACGTACAAAGTGGGGATGG;C1:37M;C2:37M375N36M;nM:0/1	255	+"
echo $myline | awk $awk_command_singlemappers

# previous test case
MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_E99_Aligned.out.bam
MW-TS-S1A-Hybr-Bst_HTMH2BGXF_S25_cat_TS.nonRibo_E99_Aligned.out.bam

# something already appears to go wrong at the intersect stage
sed -n 2082322p MW-TS-S1A-Hybr-Bst_HTMH2BGXF_S25_cat_TS.nonRibo_E99_Aligned.out.f2.singlemappers.bed

# might just be an issue because my disk was full...

awk_command_singlemappers='BEGIN {OFS="\t"; w="T"} {
    # annotate reads to know fall inbetween inton/exon boundaries, etc
    
    # jS = junction
    # 5 & 3 --> which end 
    # IN:intron
    # Cigar has splicing information
    # w: write to file (True/False)
  
    # give fields names
    chr=$1; readstart=$2; readend=$3; readname=$4; readstrand=$6; refstart=$8; refend=$9; refstrand=$10; refname=$11; genelen=$12; genestart=$13; geneend=$14;

    # obtain whether R1 or R2
    # (mw) first obtain the /1 or /2 part (if there), which becomes part of the bed files for paired reads (also remove that part)
    if (readname ~ /\/1$|\/2$/) { 
      sidx=match(readname, /\/1$|\/2$/); XM=substr(readname,sidx+1,1); readname=substr(readname,0,sidx-1) 
    } else { XM = 2 } # when theres no annotation for R1/R2, it is assumed were dealing with Read 2.
    
    # separate readname into two parts (leaving nM/NM into part 1, CG and other fields into part 2)
    sx=match(readname, /;CG:|;C1:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))
   
    if ((readstrand==refstrand && XM=2)||(readstrand!=refstrand && XM=1)) {
   
        # then check whether read fall entirely into annotation
        # if so, add "IN" to signal that
        if ((readstart >= refstart) && (readend <= refend)) {
            readname=readname";jS:IN"; w="T"
        # conversely, if both ends fall outside ref, annotate that, but in fact dont even write it (ie filter it out)
        } else if ((readstart < refstart) && (readend > refend)) {
                readname=readname";jS:OUT"; w="F";
        # if one side falls outside, annotate which part falls out
        } else if ( ((readstart < refstart)&&(readend <= refend)) || ((readstart <= refstart)&&(readend < refend)) ) {
                    if (readstrand="+") {readname=readname";jS:5"} else {readname=readname";jS:3"}; w="T"
                } else if ( ((readstart > refstart) && (readend >= refend)) || ((readstart >= refstart) && (readend > refend)) ) {
                    if (readstrand="+") {readname=readname";jS:3"} else {readname=readname";jS:5"}; w="T"
                } else {print $0 > "checkme.txt"
        }
    
        # calculate x parameter (relative position?)
        if (readstrand=="+") {
          x=1-(geneend-readend)/genelen
        } else {
          x=1-(readstart-genestart)/genelen
        }
        
        # again separate the readname into two parts
        sx=match(readname, /;CG:|;C1:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))";XM:"XM
        
        if (w=="T") {
          print chr, readstart, readend, rn, readstrand, refname, rq, refend-refstart, x}
        }
    }'