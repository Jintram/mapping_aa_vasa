Hi Martijn!

Thanks for the reminder. 

I send you the tar.gz folder with a bunch of scripts that I use to map. There are more scripts there than the ones you need, because I keep old stuff there and I am a bit lazy to clean. 

To start, I recommend you to explore the script "submit_array_starmap.sh". This is a bash script that submits different jobs in an array to perform the mapping. 
As an input, this script needs 3 arguments. The first is the name of the library that you want to map. Usually, you get 8 fastq files for each library, whose names read as:
library_L001_R1_001.fastq.gz
library_L001_R2_001.fastq.gz
library_L002_R1_001.fastq.gz
library_L002_R2_001.fastq.gz
library_L003_R1_001.fastq.gz
library_L003_R2_001.fastq.gz
library_L004_R1_001.fastq.gz
library_L004_R2_001.fastq.gz
Therefore, this first argument for the script should be for instance "library_L00", which is the name that is common in all your fastq files. 
The second argument is the protocol, usually celseq2. 
Finally, you should give the folder where you have your reference genome indexed to be mapped with STAR. 

If you then follow the different scripts that "submit_array_starmap.sh" is calling, you can get an idea of the different steps that I am following to perform the mapping. 

Good luck, and do not hesitate to ask anything you need!

All the best, 
Anna


From: Martijn Wehrens <m.wehrens@hubrecht.eu>
Date: Tuesday 24 September 2019 at 16:03
To: Anna Alemany <a.alemany@hubrecht.eu>
Subject: scripts

Hi Anna,

Thanks again for the chat!
 
Also, just as a friendly reminder, could you perhaps send me the scripts for the mapping? Many thanks in advance!
 
Best,
Martijn
 
__
Martijn Wehrens, PhD
Postdoctoral researcher
Van Rooij group, Hubrecht Institute
 
m.wehrens@hubrecht.eu
+31302121960
 