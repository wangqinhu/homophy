all:
	mkdir -p data/$(LAB)_conf
	mkdir -p data/$(LAB)_family
	./prepare_conf.pl organism.txt 2 data/seq data/$(LAB)_conf/$(GENE).conf data/$(LAB)/$(GENE).fa 1e-20
	./homophy.pl data/$(LAB)_conf/$(GENE).conf data/$(LAB)_family/$(GENE) $(DB) $(HO)
clean:
	rm -rf data/*conf/*
	rm -rf data/*family/*
	rm -rf homophy.o*
