# import numpy as np
# import pandas as pd

# filename = "lstm/dataset/all_fea_cuboctahedron_augmented.csv"

# U = np.array(
#     (
#         ((1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3)),
#         ((4, 5, 6), (4, 5, 6), (4, 5, 6), (4, 5, 6), (4, 5, 6), (4, 5, 6)),
#     )
# )
# # V = np.array(
# #     ((1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), (2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2))
# # )
# V = np.array(
#     (
#         (1, 2, 3),
#         (1, 2, 3),
#         (1, 2, 3),
#         (1, 2, 3),
#         (1, 2, 3),
#         (1, 2, 3),
#         (4, 5, 6),
#         (4, 5, 6),
#         (4, 5, 6),
#         (4, 5, 6),
#         (4, 5, 6),
#         (4, 5, 6),
#     )
# )
# K = V.copy()
# K = np.reshape(K, (2, 6, 3))
# # K.reshape((2, 6, 3))
# # K.reshape((V.shape[0] // 6, 6, V.shape[1]))
# print(K.shape)
# # print(K)
# np.random.shuffle(K)
# print(K)
# # print( // 6)
# # print(V.shape[0] % 6)
# # U = np.array((1, 2, 3, 3, 4, 5, 6, 7, 8, 9)).reshape((3, 3))
# # U.reshape((2, 6))
# # print(U.shape)
# # print(U)
# # np.random.shuffle(U)
# df = pd.read_csv(filename)
# df.to_numpy()
# print((df.shape[0]) // 101)
# print(df.shape)
# df = np.reshape(df, ((df.shape[0]) // 101, 101, df.shape[1]))
# np.random.shuffle(df)
# print(df.shape)
# # print(df.shape)
# # print(df)
# print(201 // (100 / 70), 201 // (100 / 15))
# a = df.shape[0]
# b = int(df.shape[0] // (100 / 70))
# c = int(df.shape[0] // (100 / 15)) + b + 1
# # test_len = df.shape[0] - (train_len + valid_len)
# train, validation, test = (
#     df[: b + 1, :, :],
#     df[b + 1 : c + 1, :, :],
#     df[c + 1 :, :, :],
# )
# print(train.shape, validation.shape, test.shape)


# def get_repartition_dataset(filename, train_repartition=70, validation_repartition=15):
#     df = pd.read_csv(filename)
#     df.to_numpy()
#     df = np.reshape(df, ((df.shape[0]) // 101, 101, df.shape[1]))
#     np.random.shuffle(df)

#     i = int(df.shape[0] // (100 / train_repartition))
#     j = int(df.shape[0] // (100 / validation_repartition)) + i + 1
#     train_set, validation_set, test_set = (
#         df[: i + 1, :, :],
#         df[i + 1 : j + 1, :, :],
#         df[j + 1 :, :, :],
#     )
#     return train, validation, test


from lstm.tools_database import write_shuffled_repartition_dataset
import numpy as np

# train_set, validation_set, test_set = get_repartition_dataset(
#     "lstm/dataset/all_fea_cuboctahedron_augmented.csv"
# )
# print(train_set.shape)
write_shuffled_repartition_dataset("lstm/dataset/all_fea_cuboctahedron.csv")
# U = np.array(
#     (
#         ((1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3)),
#         ((4, 5, 6), (4, 5, 6), (4, 5, 6), (4, 5, 6), (4, 5, 6), (4, 5, 6)),
#     )
# )
# V = U[:, :, :2]

# # V = V.reshape((V.shape[0], V.shape[1])
# print(V.transpose(2, 0, 1).reshape(2, 12))
