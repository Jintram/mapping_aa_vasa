#!/bin/bash

awk_command_singlemappers='BEGIN {OFS="\t"; w="T"} {
    # annotate reads to know fall inbetween inton/exon boundaries, etc
    # jS = junction
    # 5 & 3 --> which end 
    # IN:intron
    # Cigar has splicing information
    # w: write to file (True/False)
  
    # give fields names
    chr=$1; readstart=$2; readend=$3; readname=$4; readstrand=$6; refstart=$8; refend=$9; refstrand=$10; refname=$11; genelen=$12; genestart=$13; geneend=$14;
    # extract params from readname
    sx=match(readname, /;CG:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))
    
    # only look at feature when annotation is on same strand as mapped read
    if (readstrand==refstrand) {
      
       
        if ((readstart >= refstart) && (readend <= refend)) {
            readname=readname";jS:IN"; w="T"
        } else if ((readstart < refstart) && (readend > refend)) {
                readname=readname";jS:OUT"; w="F";
        } else if ( ((readstart < refstart)&&(readend <= refend)) || ((readstart <= refstart)&&(readend < refend)) ) {
                    if (readstrand="+") {readname=readname";jS:5"} else {readname=readname";jS:3"}; w="T"
                } else if ( ((readstart > refstart) && (readend >= refend)) || ((readstart >= refstart) && (readend > refend)) ) {
                    if (readstrand="+") {readname=readname";jS:3"} else {readname=readname";jS:5"}; w="T"
                } else {print $0 > "checkme.txt"
        }
    
        if (readstrand=="+") {
          x=1-(geneend-readend)/genelen
        } else {
          x=1-(readstart-genestart)/genelen
        }
        
        sx=match(readname, /;CG:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))
        
        if (w=="T") {
          print chr, readstart, readend, rn, readstrand, refname, rq, refend-refstart, x}
        }
    }'
    
    
awk_command_singlemappers_nonstranded='BEGIN {OFS="\t"; w="T"} {
    chr=$1; readstart=$2; readend=$3; readname=$4; readstrand=$6; refstart=$8; refend=$9; refstrand=$10; refname=$11; genelen=$12; genestart=$13; geneend=$14;
    sx=match(readname, /;CG:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))
    if ((readstart >= refstart) && (readend <= refend)) {
        readname=readname";jS:IN"; w="T"
    } else if ((readstart < refstart) && (readend > refend)) {
        readname=readname";jS:OUT"; w="F";
    } else if ( ((readstart < refstart)&&(readend <= refend)) || ((readstart <= refstart)&&(readend < refend)) ) {
        if (readstrand="+") {readname=readname";jS:5"} else {readname=readname";jS:3"}; w="T"
        } else if ( ((readstart > refstart) && (readend >= refend)) || ((readstart >= refstart) && (readend > refend)) ) {
        if (readstrand="+") {readname=readname";jS:3"} else {readname=readname";jS:5"}; w="T"
    } else {print $0 > "checkme.txt"}
    if (readstrand=="+") {x=1-(geneend-readend)/genelen} else {x=1-(readstart-genestart)/genelen}
    sx=match(readname, /;CG:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))
    if (w=="T") {print chr, readstart, readend, rn, readstrand, refname, rq, refend-refstart, x}
}'