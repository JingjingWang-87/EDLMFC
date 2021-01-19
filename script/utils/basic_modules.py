from keras.layers import Dense, Conv1D, concatenate, BatchNormalization, MaxPooling1D, Bidirectional, LSTM
from keras.layers import Dropout, Flatten, Input
from keras.models import Model

def conjoint_struct_cnn_blstm(pro_coding_length, rna_coding_length, vector_repeatition_cnn):

    if type(vector_repeatition_cnn)==int:
        vec_len_p = vector_repeatition_cnn
        vec_len_r = vector_repeatition_cnn
    else:
        vec_len_p = vector_repeatition_cnn[0]
        vec_len_r = vector_repeatition_cnn[1]

    # NN for protein feature analysis by one hot encoding
    xp_in_conjoint_struct_cnn_blstm = Input(shape=(pro_coding_length, vec_len_p))
    xp_cnn = Conv1D(filters=45, kernel_size=6, strides=1, activation='relu')(xp_in_conjoint_struct_cnn_blstm)
    xp_cnn = MaxPooling1D(pool_size=2)(xp_cnn)
    xp_cnn = BatchNormalization()(xp_cnn)
    xp_cnn = Dropout(0.2)(xp_cnn)
    xp_cnn = Conv1D(filters=64, kernel_size=6, strides=1, activation='relu')(xp_cnn)
    xp_cnn = MaxPooling1D(pool_size=2)(xp_cnn)
    xp_cnn = BatchNormalization()(xp_cnn)
    xp_cnn = Dropout(0.2)(xp_cnn)
    xp_cnn = Conv1D(filters=86, kernel_size=6, strides=1, activation='relu')(xp_cnn)
    xp_cnn = BatchNormalization()(xp_cnn)
    xp_cnn = Dropout(0.2)(xp_cnn)
    xp_cnn = Bidirectional(LSTM(45,return_sequences=True))(xp_cnn)
    xp_cnn = Flatten()(xp_cnn)
    xp_out_conjoint_cnn_blstm = Dense(64)(xp_cnn)
    xp_out_conjoint_cnn_blstm = Dropout(0.2)(xp_out_conjoint_cnn_blstm)

    # NN for RNA feature analysis  by one hot encoding
    xr_in_conjoint_struct_cnn_blstm = Input(shape=(rna_coding_length, vec_len_r))
    xr_cnn = Conv1D(filters=45, kernel_size=6, strides=1, activation='relu')(xr_in_conjoint_struct_cnn_blstm)
    xr_cnn = MaxPooling1D(pool_size=2)(xr_cnn)
    xr_cnn = BatchNormalization()(xr_cnn)
    xr_cnn = Dropout(0.2)(xr_cnn)
    xr_cnn = Conv1D(filters=64, kernel_size=5, strides=1, activation='relu')(xr_cnn)
    xr_cnn = MaxPooling1D(pool_size=2)(xr_cnn)
    xr_cnn = BatchNormalization()(xr_cnn)
    xr_cnn = Dropout(0.2)(xr_cnn)
    xr_cnn = Conv1D(filters=86, kernel_size=5, strides=1, activation='relu')(xr_cnn)
    xr_cnn = MaxPooling1D(pool_size=2)(xr_cnn)
    xr_cnn = BatchNormalization()(xr_cnn)
    xr_cnn = Dropout(0.2)(xr_cnn)
    xr_cnn = Bidirectional(LSTM(45,return_sequences=True))(xr_cnn)
    xr_cnn = Flatten()(xr_cnn)
    xr_out_conjoint_cnn_blstm = Dense(64)(xr_cnn)
    xr_out_conjoint_cnn_blstm = Dropout(0.2)(xr_out_conjoint_cnn_blstm)


    x_out_conjoint_cnn_blstm = concatenate([xp_out_conjoint_cnn_blstm, xr_out_conjoint_cnn_blstm])
    x_out_conjoint_cnn_blstm = Dense(128, kernel_initializer='random_uniform', activation='relu')(x_out_conjoint_cnn_blstm)
    x_out_conjoint_cnn_blstm = Dropout(0.25)(x_out_conjoint_cnn_blstm)
    x_out_conjoint_cnn_blstm = BatchNormalization()(x_out_conjoint_cnn_blstm)
    x_out_conjoint_cnn_blstm = Dropout(0.3)(x_out_conjoint_cnn_blstm)
    x_out_conjoint_cnn_blstm = Dense(64, kernel_initializer='random_uniform', activation='relu')(x_out_conjoint_cnn_blstm)
    y_conjoint_cnn_blstm = Dense(2, activation='softmax')(x_out_conjoint_cnn_blstm)

    model_conjoint_struct_cnn_blstm = Model(inputs=[xp_in_conjoint_struct_cnn_blstm, xr_in_conjoint_struct_cnn_blstm], outputs=y_conjoint_cnn_blstm)


    return model_conjoint_struct_cnn_blstm
