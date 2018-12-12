FROM pshannon/trenavizbase
USER root
COPY Rprofile.site /usr/local/lib/R/etc/
USER trena
WORKDIR /home/trena/app
COPY packages/*gz ./
COPY igap.R ./
COPY testDB.R ./
RUN R CMD INSTALL --no-test-load -l ~/library trena_1.5.10.tar.gz
RUN R CMD INSTALL --no-test-load -l ~/library TrenaProject_0.99.24.tar.gz
RUN R CMD INSTALL --no-test-load -l ~/library TrenaProjectIGAP_0.99.15.tar.gz
RUN R CMD INSTALL --no-test-load -l ~/library TrenaViz_0.99.10.tar.gz
EXPOSE 3838
CMD ["R", "-f", "igap.R"]





