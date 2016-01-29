all:
	./prepare_conf.pl organism.txt 2 data/seq data/conf/$(GENE).conf data/seed/$(GENE).fa 1e-20
	./homophy.pl data/conf/$(GENE).conf data/family/$(GENE)
clean:
	rm -rf data/conf/*
	rm -rf data/family/*
	rm -rf homophy.o*
