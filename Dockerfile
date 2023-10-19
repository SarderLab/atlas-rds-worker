FROM satijalab/seurat:latest

RUN R -e "install.packages('argparser', dependencies=TRUE, repos='http://cran.rstudio.com/')"

# RUN R -e "remotes::install_github('satijalab/seurat', 'seurat5', quiet = TRUE)"
# RUN R -e "remotes::install_github('satijalab/seurat-data', 'seurat5', quiet = TRUE)"
# RUN R -e "remotes::install_github('satijalab/azimuth', 'seurat5', quiet = TRUE)"
# RUN R -e "remotes::install_github('satijalab/seurat-wrappers', 'seurat5', quiet = TRUE)"
# RUN R -e "remotes::install_github('stuart-lab/signac', 'seurat5', quiet = TRUE)"
# RUN R -e "remotes::install_github('mojaveazure/seurat-object', 'seurat5', quiet = TRUE)"

RUN mkdir /project
WORKDIR /project

# RUN wget https://data.kitware.com/api/v1/file/hashsum/sha512/f9242aae27bfe2d1ce190a83c322abe5f1ee2e830bba5012f136a98ab3cec219825cc616b6b3e4f04f6641ee5c048747d20603032681ac9272e40653e078a1cf/download -o KidneyAtlas_snCV3_20percent.h5Seurat

COPY . /project/.
WORKDIR /project/cli

RUN cat slicer_cli_list.json

ENTRYPOINT ["/bin/bash", "docker-entrypoint.sh"]
