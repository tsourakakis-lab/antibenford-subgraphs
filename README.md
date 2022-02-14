# Anomaly detection in Financial Networks 
## Tianyi Chen, Charalampos E. Tsourakakis

The experimental part is written in combination of Python3(Jupyter notebooks) and C++, to achieve easy demonstration and efficiency simultaneously. 

Practitioners can first go through the notebooks, which preprocess the raw data and postprocess the DSD outputs. The DSD component in DSD_cpp is implemented in C++, please see the README inside the folder for compile instruction.

Our datasets are available on Google drive: https://drive.google.com/drive/folders/1NKjzJS7w1dDqvwMm-lQcP3ryf3X8AfGx?usp=sharing

Two folders are included in the pointed drive, ETH-anomaly and ETH-scalability. The first folder contains raw data for ETH-Jan-2018 and ETH-Jan-2019. Users can preprocess it with the corresponding notebooks in this repository. The second folder contains five prerpocessed networks for scalability evaluation.

Finally, we omitted some less important information and take monthly range when construct our network, as the whole dataset is too huge to dump. We encourage anyone interested to varify our findings and explore more on https://www.kaggle.com/bigquery/ethereum-blockchain.
