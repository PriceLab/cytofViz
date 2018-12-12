
ID = 08dec2018
build:
	docker build -t trenaigap:$(ID) -f Dockerfile .

bash:
	docker run -ti --rm -p 60051:3838 trenaigap:$(ID) bash

dockerRun:
	docker run -p 60051:3838 trenaigap:$(ID)

dockerView:
	open http://localhost:60051

save:
	time docker save trenaigap:$(ID) >  ~/tmp/dockerIGAP.$(ID)
	ls -lat ~/tmp/dockerIGAP.$(ID)

scp:
	time scp ~/tmp/dockerIGAP.$(ID) pshannon@trena:/local/tmp/
	ssh pshannon@trena ls -l /local/tmp/dockerIGAP.$(ID)
#	ssh pshannon@trena docker load < /local/tmp/dockerIGAP.$(ID)


saveTest:
	time docker save trenaigap:test >  ~/tmp/dockerIGAP.test.$(ID)
	ls -lat ~/tmp/dockerIGAP.test.$(ID)


scpTest:
	time scp ~/tmp/dockerIGAP.test.$(ID) pshannon@trena:/local/tmp/
	ssh pshannon@trena ls -l /local/tmp/dockerIGAP.test.$(ID)


cleanDocker:
	  docker rm $$(docker ps -q -f 'status=exited')
	  docker rmi $$(docker images -q -f "dangling=true")
