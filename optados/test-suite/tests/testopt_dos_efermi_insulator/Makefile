all: $(patsubst %.bands.bz2, %.bands, $(wildcard *.bands.bz2))

%.bands: %.bands.bz2
	$(eval SEEDNAME:=$(patsubst %.bands.bz2, %, $<))
	cat $< | bunzip2 - > $(SEEDNAME).bands

.PHONY: all
