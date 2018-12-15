ID = 27nov2018
build:
	docker build -t cytof -f Dockerfile .

bash:
	docker run -ti --rm -p 60051:3838 cytof bash

dockerRun:
	docker run -p 60051:3838 cytof

save:
	time docker save cytof >  ~/tmp/docker_cytof.image
	ls -lat ~/tmp/docker_cytof.image

scp:
	time scp ~/tmp/docker_cytof.image pshannon@trena:/local/tmp/
	ssh pshannon@trena ls -l /local/tmp/docker_cytof.image
#	ssh pshannon@trena docker load < /local/tmp/dockerIGAP.$(ID)


saveTest:
	time docker save trenaigap:test >  ~/tmp/dockerIGAP.test.$(ID)
	ls -lat ~/tmp/dockerIGAP.test.$(ID)


scpTest:
	time scp ~/tmp/dockerIGAP.test.$(ID) pshannon@trena:/local/tmp/
	ssh pshannon@trena ls -l /local/tmp/dockerIGAP.test.$(ID)

