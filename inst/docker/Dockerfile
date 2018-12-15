FROM pshannon/trenashinybase
USER root
USER trena
WORKDIR /home/trena/app
COPY cytofViz.R .
COPY tbl.allDays.allProteins.RData .
CMD ["R", "-e", "shiny::runApp('/home/trena/app/cytofViz.R', port=3838, host='0.0.0.0')"]





