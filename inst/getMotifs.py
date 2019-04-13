#!/usr/bin/python
# -*- coding: utf-8 -*-
"""636
Created on Thu Apr  4 20:41:25 2019

@author: Jun Wang
"""
def DeepCirCode_getMotifs(x_train,y_train, x_test, y_test):
    import numpy as np

    x_train = np.array(x_train)
    y_train = np.array(y_train)
    x_test = np.array(x_test)
    y_test = np.array(y_test)

    import keras
    from keras import backend as K
    from keras.models import Sequential
    from keras.layers import Conv1D, MaxPooling1D
    from keras.layers import Dense, Dropout, Flatten
    from keras import initializers
    from keras import optimizers
    from keras import losses

    np.random.seed(123)

    model = Sequential()

    conv_layer1 = Conv1D(filters=128,
                    kernel_size=12,
                    strides=1,
                    padding='valid',
                    activation='relu',
                    input_shape=(200,4))

    conv_layer2 = Conv1D(filters=128,
                    kernel_size=6,
                    strides=1,
                    padding='valid',
                    activation='relu')

    model.add(conv_layer1)
    model.add(Dropout(0))

    model.add(conv_layer2)
    model.add(Dropout(0.5))

    model.add(MaxPooling1D(pool_size=4, strides=4))
    model.add(Dropout(0.5))

    model.add(Flatten())
    model.add(Dense(2,activation='softmax'))

    model.summary()

    model.compile(loss=losses.binary_crossentropy,
              optimizer=optimizers.rmsprop(),
              metrics=['accuracy'])

    model.fit(x_train, y_train,
          batch_size=128,
          epochs=80,
          validation_split=0.1)

    # motif visualization:
    ### motif visualization for BS model human:
    conv_output = conv_layer1.get_output_at(0)
    f = K.function([model.input], [K.argmax(conv_output, axis=1), K.max(conv_output, axis=1)])

    motifs = np.zeros((128, 12, 4))
    nsites = np.zeros(128)

    # select the positive samples from test set:
    row = y_test[:,1]==1.0
    y_test_positive = y_test[np.ix_(row)]
    x_test_positive = x_test[np.ix_(row)]

    for i in range(0, len(x_test_positive), 100):
        x = x_test_positive[i:i+100]
        z = f([x])
        max_inds = z[0] # N x M matrix, where M is the number of motifs
        max_acts = z[1]
        for m in range(128):
            for n in range(len(x)):
                # Forward strand
                if max_acts[n, m] > 0:
                    nsites[m] += 1
                    motifs[m] += x[n, max_inds[n, m]:max_inds[n, m] + 12, :]

    print('Making motifs')
    motifs = motifs[:, :, [0, 3, 2, 1]]
    motifs_file = open('DeepCirCode_Position_Probability_Matrix.txt', 'w')
    motifs_file.write('MEME version 4.9.0\n\n'
                  'ALPHABET= ACGU\n\n'
                  'strands: + -\n\n'
                  'Background letter frequencies (from uniform background):\n'
                  'A 0.25000 C 0.25000 G 0.25000 U 0.25000\n\n')

    for m in range(128):
        if nsites[m] == 0:
            continue
        motifs_file.write('MOTIF M_n%i O%i\n' % (m, m))
        motifs_file.write("letter-probability matrix: alength= 4 w= %i nsites= %i E= 1337.0e-6\n" % (12, nsites[m]))
        for j in range(12):
            motifs_file.write("%f %f %f %f\n" % tuple(1.0 * motifs[m, j, 0:4] / np.sum(motifs[m, j, 0:4])))
        motifs_file.write('\n')

    motifs_file.close()
    print("Done")

