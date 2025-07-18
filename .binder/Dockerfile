FROM rocker/binder:latest
LABEL maintainer='Ka Yin Leung'
COPY --chown=${NB_USER} . ${HOME}
USER ${NB_USER}



RUN wget https://github.com/rivm-syso/mpox-network/raw/main/DESCRIPTION \
&& R -e "devtools::install_deps(dependencies = TRUE)" \
&& R -e "install.packages('pak')" \
&& R -e "pak::pak('rivm-syso/mpox-network')" 

RUN rm DESCRIPTION.1; exit 0
