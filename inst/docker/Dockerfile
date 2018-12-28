FROM pshannon/trenashinybase
USER root
WORKDIR /home/trena/app
RUN mkdir /home/trena/app/images
COPY cytofViz.R .
COPY tbl.allDays.allProteins.transformed.RData .
COPY images/phenograph.png images/
USER trena
CMD ["R", "-e", "shiny::runApp('/home/trena/app/cytofViz.R', port=3838, host='0.0.0.0')"]





