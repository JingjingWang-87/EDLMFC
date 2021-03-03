# EDLMFC
 An ensemble deep learning with multi-scale feature  combination for ncRNA-protein interaction prediction. 

The _untils_, _data_ and _result_ directories contain model codes, tested data sets and generated results, respectively.
The depended python packages are listed in _requirements.txt_. The package versions should be followed by users in their environments to achieve the supposed performance.

## How to run

The program is in Python 3.6 using [Keras](https://keras.io/) and [Tensorflow](https://www.tensorflow.org/) backends. Use the below bash command to run EDLMFC.

```bash
    python main.py -d dataset
```

The parameter of _dataset_ could be RPI488, RPI1807 and NPInter v2.0. Then, EDLMFC will perform 5-fold cross validation on the specific dataset.


## Two RPI datasets

The widely used RPI benchmark datasets are organized in the _data_ directory. 

Due to the limitation of the hardware conditions of the selected RNA secondary structure method, it can only predict the secondary structure of RNA with a length of no more than 500 nucleotides, so we preprocessed the data.

                       Dataset    | #Positive pairs | #Negative pairs | RNAs | Proteins |Reference
-------------|--------------|-----------------|------------------|--------|----------|------------

Original set     RPI1807             1807                  1436            1078       3131        [1]
                 NPInter v2.0        10412                10412           4636        449         [2]             
                 RPI488              243                   245             25         247         [1]
-------------|--------------|-----------------|------------------|--------|----------|-------------

Optimal set     RPI1807              652                    221              646        868         [1]
                NPInter v2.0|        1943                   1943            513        448          [2]       
                RPI488                43                    233             13          155         [1]
## Help

For any questions, feel free to contact me by WangJingJing@emails.bjut.edu.cn or start an issue instead.


[1] Pan, X.Y.; Fan, Y.X.; Yan, J.C.; Shen, H.B. IPMiner: hidden ncRNA-protein interaction sequential pattern mining with stacked autoencoder for accurate computational prediction. Bmc Genomics 2016, 17. doi:ARTN 582 10.1186/s12864-016-2931-8.

[2] Yuan, J.;Wu,W.; Xie, C.Y.; Zhao, G.G.; Zhao, Y.; Chen, R.S. NPInter v2.0: an updated database of ncRNA interactions. Nucleic Acids Research 2014, 42, D104–D108. doi:10.1093/nar/gkt1057.
