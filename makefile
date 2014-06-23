CC=g++
LIB=libabre.a

CINCLUDE=-I../trollbase
#COPT=-g -Wno-deprecated
COPT=-ffast-math -O3 -DUNIXVER -Wno-deprecated -Wall
#COPT=-march=i686 -g -DUNIXVER -Wno-deprecated -static -Wall
CLIB=
#CLIB=../trollbase/libTroll.a

# removed struct frag* and gn2* object files

$(LIB): 	aa_seq.o alib.o ali_dist.o ali_strand_eval.o ali_str_info.o application.o \
		argv.o clusterset.o dpmatrix.o fastaio.o formats.o \
		gstrings.o hmap_eval.o hmap2_eval.o hmapalib_seq.o \
		hmapio.o kmedoidclusterer.o \
		noalib.o pirio.o pstore.o rcfile.o sequence.o sflags.o skel_ali.o \
		skel_set.o submatrix.o UPGMA_Clusterer.o \
		UPGMA_Tree.o
	    	ar rv $(LIB) $?

aa_seq.o:	aa_seq.cpp aa_seq.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

#akalib.o:	akalib.cpp akalib.h
#		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

alib.o:		alib.cpp alib.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

ali_dist.o:	ali_dist.cpp ali_dist.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

ali_frag.o:	ali_frag.cpp ali_frag.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

ali_strand_eval.o:	ali_strand_eval.cpp ali_strand_eval.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

ali_str_info.o:	ali_str_info.cpp ali_str_info.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

application.o:	application.cpp application.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

argv.o:		argv.cpp argv.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

dpmatrix.o:	dpmatrix.cpp dpmatrix.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

fastaio.o:	fastaio.cpp fastaio.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

formats.o:	formats.cpp formats.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

frag_matrix.o:	frag_matrix.cpp frag_matrix.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

frag_set.o:	frag_set.cpp frag_set.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

gn2_eval.o:	gn2_eval.cpp gn2_eval.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

gnoalib.o:	gnoalib.cpp gnoalib.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

gn2lib_seq.o:	gn2lib_seq.cpp gn2lib_seq.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

#gnoalib_seq.o:	gnoalib_seq.cpp gnoalib_seq.h
#		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

gstrings.o:	gstrings.cpp gstrings.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

hmap_eval.o:	hmap_eval.cpp hmap_eval.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

hmap2_eval.o:	hmap2_eval.cpp hmap_eval.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

hmapalib_seq.o:	hmapalib_seq.cpp hmapalib_seq.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

hmapio.o:	hmapio.cpp hmapio.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

kmedoidclusterer.o:	kmedoidclusterer.cpp kmedoidclusterer.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

noalib.o:	noalib.cpp noalib.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

pirio.o:	pirio.cpp pirio.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

pstore.o:	pstore.cpp pstore.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

rcfile.o:	rcfile.cpp rcfile.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

sequence.o: 	sequence.cpp sequence.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

sflags.o:	sflags.cpp sflags.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

skel_ali.o:	skel_ali.cpp skel_ali.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

skel_set.o:	skel_set.cpp skel_set.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

clusterset.o:	clusterset.cpp clusterset.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

sse_frag_set.o:	sse_frag_set.cpp sse_frag_set.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

struct.o:	struct.cpp struct.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

submatrix.o:	submatrix.cpp submatrix.h
		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

#UPGMA_Clusterer.o:	UPGMA_Clusterer.cpp UPGMA_Clusterer.h
#		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

#UPGMA_Tree.o:	UPGMA_Tree.cpp UPGMA_Tree.h
#		$(CC) -c $(COPT) -o $@ $(CINCLUDE) $(@F:.o=.cpp)

aaa:		aa_ali.cpp $(LIB)
		$(CC) $(COPT) -o $@ $(CINCLUDE) aa_ali.cpp $(LIB) $(CLIB)

cn_acc_analys:	cn_acc_analys.cpp $(LIB)
		$(CC) $(COPT) -o $@ $(CINCLUDE) $(@F:=.cpp) $(LIB) $(CLIB)

gn2:		gn2.cpp gn2_eval.h cw.h ucw.h kscw.h crcw.h $(LIB)
		$(CC) $(COPT) -o $@ $(CINCLUDE) $(@F:=.cpp) $(LIB) $(CLIB)

gnoali:		gnoali.cpp $(LIB)
		$(CC) $(COPT) -o $@ $(CINCLUDE) $(@F:=.cpp) $(LIB) $(CLIB)

Frag_Analyzer:	Frag_Analyzer.cpp $(LIB)
		$(CC) $(COPT) -o $@ $(CINCLUDE) $(@F:=.cpp) $(LIB) $(CLIB)

S4_align:	S4_align.cpp $(LIB)
		$(CC) $(COPT) -o $@ $(CINCLUDE) $(@F:=.cpp) $(LIB) $(CLIB)

S4_align_gn2:	S4_align_gn2.cpp $(LIB)
		$(CC) $(COPT) -o $@ $(CINCLUDE) $(@F:=.cpp) $(LIB) $(CLIB)

get_shifts:	get_shifts.cpp alignment.h $(LIB)
		$(CC) $(COPT) -o $@ $(CINCLUDE) $(@F:=.cpp) $(LIB) $(CLIB)

S4_one_ali:	S4_one_ali.cpp $(LIB)
		$(CC) $(COPT) -o $@ $(CINCLUDE) $(@F:=.cpp) $(LIB) $(CLIB)

get_area_diffs:	get_area_diffs.cpp $(LIB)
		$(CC) $(COPT) -o $@ $(CINCLUDE) $(@F:=.cpp) $(LIB) $(CLIB)

get_loop_dists:	get_loop_dists.cpp $(LIB)
		$(CC) $(COPT) -o $@ $(CINCLUDE) $(@F:=.cpp) $(LIB) $(CLIB)

get_loop_lengths:	get_loop_lengths.cpp $(LIB)
		$(CC) $(COPT) -o $@ $(CINCLUDE) $(@F:=.cpp) $(LIB) $(CLIB)

nalign:		nalign.cpp hmap_eval.h cw.h ucw.h $(LIB)
		$(CC) $(COPT) -o $@ $(CINCLUDE) $(@F:=.cpp) $(LIB) $(CLIB)

nalign2:	nalign2.cpp hmap2_eval.h cw.h ucw.h kscw.h crcw.h $(LIB)
		$(CC) $(COPT) -o $@ $(CINCLUDE) $(@F:=.cpp) $(LIB) $(CLIB)

clean:          
		rm -f *.o $(LIB)

distclean:
		rm -f *.o *~ $(LIB)
