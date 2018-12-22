ID = 21dec2018
build:
	docker build -t cytof_$(ID) -f Dockerfile .

bash:
	docker run -ti --rm -p 60051:3838 cytof_$(ID) bash

dockerRun:
	docker run -p 60051:3838 cytof_$(ID)

save:
	time docker save cytof_$(ID) >  ~/tmp/docker_cytof_$(ID).image
	ls -lat ~/tmp/docker_cytof_$(ID).image

scp:
	time scp ~/tmp/docker_cytof_$(ID).image pshannon@trena:/local/tmp/
	ssh pshannon@trena ls -l /local/tmp/docker_cytof_$(ID).image
#	ssh pshannon@trena docker load < /local/tmp/docker_cytof_$(ID).image


