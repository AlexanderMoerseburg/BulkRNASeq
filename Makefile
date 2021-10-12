REFERENCES= /rds/project/rds-O11U8YqSuCk/SM/annotations 

Homo_sapiens_Ensembl_GRCh37.tar.gz:
	cd ${REFERENCES}
	wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz
	tar -xzvf GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz
	ls Homo_sapiens/Ensembl/GRCh37/Annotation 
	ls Homo_sapiens/Ensembl/GRCh37/Sequence

Mus_musculus_Ensembl_GRCm38.tar.gz:
	cd ${REFERENCES}
	wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/Ensembl/GRCm38/Mus_musculus_Ensembl_GRCm38.tar.gz
	tar -xzvf Mus_musculus_Ensembl_GRCm38.tar.gz
	ls Mus_musculus/Ensembl/GRCm38/Annotation
	ls Mus_musculus/Ensembl/GRCm38/Sequence 
