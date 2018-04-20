## Image processing and primary analysis on germination process
1. You need to download loci_tools and 'imreadBF' and 'imreadBFmeta' to parse the nd2 file. 
2. Start new matlab script, define workpath of nd2 movie file and call funciton 'score_all11.m' (which inturn call the other embedded functions 'score_well_red11' and 'score_field_of_view_red11'. This will identify individual spores (or cluster if not correctly separated) and track each germination process by generating movie strips (frames of images) for each spore/cluster and caculate features (e.g. spore size, major axis length) to give rough estimation on germination time. Movie strips are saved as 'IMGs11.m' files while primary analysis are saved as 'results11.m'.

3. Run 'MANUAL_scoring_doublet11' to mannualy curate the primary analysis from previous steps by 
(1) identifying spores touching in pairs (doublet) and calling function 'movie_strip_score_doublet_red_MANUALLY' to bring up interactive popup windows that allows clicking on the right frame of germination. 
(2) identifying isolated spores (singles) and calling function 'movie_strip_score_germination_red_MANUALLY' to bring up interactive popup windows that allows clicking on the right frame of germination. 
Analysis results are added to the previous results and saved as 'results_with_manual_doublets'