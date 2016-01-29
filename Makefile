all:
	./prepare_conf.pl organism.txt 2 data/seq data/conf/$(GENE).conf data/seed/$(GENE).fa
	./homophy.pl data/conf/$(GENE).conf $(GENE)
clean:
	rm -rf $(GENE)
	rm formatdb.log
