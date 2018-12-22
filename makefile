all:  docs install

docs:
	R -e "devtools::document()"

build:
	(cd ..; R CMD build --no-build-vignettes trenaViz)

install:
	(cd ..; R CMD INSTALL --no-test-load trenaViz)

test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done

