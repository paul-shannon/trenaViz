build: webpack assemble buildPackage installPackage

assemble:
	(cd ./dist; m4 -P trenaviz.html-template | sed s/\`//g > trenaviz.html)
	ls -l ./dist/trenaviz.html

buildPackage:
	(cd ../../..; R CMD build --no-build-vignettes trenaViz)

installPackage:
	(cd ../../..; R CMD INSTALL `ls -at trenaViz_* | head -1`)

webpack:
	npm update
	webpack

test:
	(cd ./dist; pwd)

all: clean install run

clean:
	- rm -rf dist/
	- rm -rf node_modules

install:
	npm install
