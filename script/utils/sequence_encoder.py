
import numpy as np

# encoder for protein sequence
class ProEncoder:
    elements = 'AIYHRDC'
    structs = 'HEC'

    element_number = 7
    # number of structure kind
    struct_kind = 3

    # clusters: {A,G,V}, {I,L,F,P}, {Y,M,T,S}, {H,N,Q,W}, {R,K}, {D,E}, {C}
    pro_intab = 'AGVILFPYMTSHNQWRKDEC'
    pro_outtab = 'AAAIIIIYYYYHHHHRRDDC'

    def __init__(self, WINDOW_P_UPLIMIT, WINDOW_P_STRUCT_UPLIMIT, CODING_FREQUENCY, VECTOR_REPETITION_CNN,
                 TRUNCATION_LEN=None, PERIOD_EXTENDED=None):

        self.WINDOW_P_UPLIMIT = WINDOW_P_UPLIMIT  ##self自动上传变量
        self.WINDOW_P_STRUCT_UPLIMIT = WINDOW_P_STRUCT_UPLIMIT
        self.CODING_FREQUENCY = CODING_FREQUENCY
        self.VECTOR_REPETITION_CNN = VECTOR_REPETITION_CNN

        self.TRUNCATION_LEN = TRUNCATION_LEN
        self.PERIOD_EXTENDED = PERIOD_EXTENDED

        # list and position map for k_mer
        k_mers = ['']
        self.k_mer_list = []
        self.k_mer_map = {}
        for T in range(self.WINDOW_P_UPLIMIT):
            temp_list = []
            for k_mer in k_mers:
                for x in self.elements:
                    temp_list.append(k_mer + x)
            k_mers = temp_list
            self.k_mer_list += temp_list
        for i in range(len(self.k_mer_list)):
            self.k_mer_map[self.k_mer_list[i]] = i

        # list and position map for k_mer structure
        k_mers = ['']
        self.k_mer_struct_list = []
        self.k_mer_struct_map = {}
        for T in range(self.WINDOW_P_STRUCT_UPLIMIT):
            temp_list = []
            for k_mer in k_mers:
                for s in self.structs:
                    temp_list.append(k_mer + s)
            k_mers = temp_list
            self.k_mer_struct_list += temp_list
        for i in range(len(self.k_mer_struct_list)):
            self.k_mer_struct_map[self.k_mer_struct_list[i]] = i

        # table for amino acid clusters
        self.transtable = str.maketrans(self.pro_intab, self.pro_outtab)


        # print(len(self.k_mer_list))
        # print(self.k_mer_list)
        # print(len(self.k_mer_struct_list))
        # print(self.k_mer_struct_list)

    def encode_conjoint(self, seq):
        seq = seq.translate(self.transtable)
        seq = ''.join([x for x in seq if x in self.elements])
        seq_len = len(seq)
        if seq_len == 0:
            return 'Error'
        result = []
        offset = 0
        for K in range(1, self.WINDOW_P_UPLIMIT + 1):
            vec = [0.0] * (self.element_number ** K)
            counter = seq_len - K + 1
            for i in range(counter):
                k_mer = seq[i:i + K]
                vec[self.k_mer_map[k_mer] - offset] += 1
            vec = np.array(vec)
            offset += vec.size
            if self.CODING_FREQUENCY:
                vec = vec / vec.max()
            result += list(vec)
        return np.array(result)

    def encode_conjoint_struct(self, seq, struct):
        seq = seq.translate(self.transtable)
        seq_temp = []
        struct_temp = []
        for i in range(len(seq)):
            if seq[i] in self.elements:
                seq_temp.append(seq[i])
                struct_temp.append(struct[i])
        seq = ''.join(seq_temp)
        struct = ''.join(struct_temp)
        seq_len = len(seq)
        if seq_len == 0:
            return 'Error'

        result_seq = []
        offset_seq = 0
        for K in range(1, self.WINDOW_P_UPLIMIT + 1):
            vec_seq = [0.0] * (self.element_number ** K)
            counter = seq_len - K + 1
            for i in range(counter):
                k_mer = seq[i:i + K]
                vec_seq[self.k_mer_map[k_mer] - offset_seq] += 1
            vec_seq = np.array(vec_seq)
            offset_seq += vec_seq.size
            if self.CODING_FREQUENCY:
                vec_seq = vec_seq / vec_seq.max()
            result_seq += list(vec_seq)


        result_struct = []
        offset_struct = 0
        for K in range(1, self.WINDOW_P_STRUCT_UPLIMIT + 1):
            vec_struct = [0.0] * (self.struct_kind ** K)
            counter = seq_len - K + 1
            for i in range(counter):
                k_mer_struct = struct[i:i + K]
                vec_struct[self.k_mer_struct_map[k_mer_struct] - offset_struct] += 1
            vec_struct = np.array(vec_struct)
            offset_struct += vec_struct.size
            if self.CODING_FREQUENCY:
                vec_struct = vec_struct / vec_struct.max()
            result_struct += list(vec_struct)
        return np.array(result_seq + result_struct)


# encoder for RNA sequence
class RNAEncoder:
    elements = 'AUCG'
    structs = 'SMHBEIX'

    element_number = 4
    struct_kind = 7

    def __init__(self, WINDOW_R_UPLIMIT, WINDOW_R_STRUCT_UPLIMIT, CODING_FREQUENCY, VECTOR_REPETITION_CNN,
                 TRUNCATION_LEN=None, PERIOD_EXTENDED=None):

        self.WINDOW_R_UPLIMIT = WINDOW_R_UPLIMIT
        self.WINDOW_R_STRUCT_UPLIMIT = WINDOW_R_STRUCT_UPLIMIT
        self.CODING_FREQUENCY = CODING_FREQUENCY
        self.VECTOR_REPETITION_CNN = VECTOR_REPETITION_CNN

        self.TRUNCATION_LEN = TRUNCATION_LEN
        self.PERIOD_EXTENDED = PERIOD_EXTENDED

        # list and position map for k_mer
        k_mers = ['']
        self.k_mer_list = []
        self.k_mer_map = {}
        for T in range(self.WINDOW_R_UPLIMIT):
            temp_list = []
            for k_mer in k_mers:
                for x in self.elements:
                    temp_list.append(k_mer + x)
            k_mers = temp_list
            self.k_mer_list += temp_list
        for i in range(len(self.k_mer_list)):
            self.k_mer_map[self.k_mer_list[i]] = i

        # list and position map for k_mer structure
        k_mers = ['']
        self.k_mer_struct_list = []
        self.k_mer_struct_map = {}
        for T in range(self.WINDOW_R_STRUCT_UPLIMIT):
            temp_list = []
            for k_mer in k_mers:
                for s in self.structs:
                    temp_list.append(k_mer + s)
            k_mers = temp_list
            self.k_mer_struct_list += temp_list
        for i in range(len(self.k_mer_struct_list)):
            self.k_mer_struct_map[self.k_mer_struct_list[i]] = i

        # print(len(self.k_mer_list))
        # print(self.k_mer_list)
        # print(len(self.k_mer_struct_list))
        # print(self.k_mer_struct_list)

    def encode_conjoint(self, seq):
        seq = seq.replace('T', 'U')
        seq = ''.join([x for x in seq if x in self.elements])
        seq_len = len(seq)
        if seq_len == 0:
            return 'Error'
        result = []
        offset = 0
        for K in range(1, self.WINDOW_R_UPLIMIT + 1):
            vec = [0.0] * (self.element_number ** K)
            counter = seq_len - K + 1
            for i in range(counter):
                k_mer = seq[i:i + K]
                vec[self.k_mer_map[k_mer] - offset] += 1
            vec = np.array(vec)
            offset += vec.size
            if self.CODING_FREQUENCY:
                vec = vec / vec.max()
            result += list(vec)
        return np.array(result)

    def encode_conjoint_struct(self, seq, struct):
        seq = seq.replace('T', 'U')
        struct = struct.replace(')', '(')

        seq_temp = []
        struct_temp = []
        for i in range(len(seq)):
            if seq[i] in self.elements:
                seq_temp.append(seq[i])
                struct_temp.append(struct[i])
        seq = ''.join(seq_temp)
        struct = ''.join(struct_temp)
        seq_len = len(seq)
        if seq_len == 0:
            return 'Error'

        result_seq = []
        offset_seq = 0
        for K in range(1, self.WINDOW_R_UPLIMIT + 1):
            vec_seq = [0.0] * (self.element_number ** K)
            counter = seq_len - K + 1
            for i in range(counter):
                k_mer = seq[i:i + K]
                vec_seq[self.k_mer_map[k_mer] - offset_seq] += 1
            vec_seq = np.array(vec_seq)
            offset_seq += vec_seq.size
            if self.CODING_FREQUENCY:
                vec_seq = vec_seq / vec_seq.max()
            result_seq += list(vec_seq)


        result_struct = []
        offset_struct = 0
        for K in range(1, self.WINDOW_R_STRUCT_UPLIMIT + 1):
            vec_struct = [0.0] * (self.struct_kind ** K)
            counter = seq_len - K + 1
            for i in range(counter):
                k_mer_struct = struct[i:i + K]
                vec_struct[self.k_mer_struct_map[k_mer_struct] - offset_struct] += 1
            vec_struct = np.array(vec_struct)
            offset_struct += vec_struct.size
            if self.CODING_FREQUENCY:
                vec_struct = vec_struct / vec_struct.max()
            result_struct += list(vec_struct)
        return np.array(result_seq + result_struct)
