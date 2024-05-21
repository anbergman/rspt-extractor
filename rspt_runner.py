import numpy as np
import rspt_exchange as rspt


# def extract_projections(exchange):
#     natom = exchange.shape[0]
#     j_dict = {}
#     j_dict["scalar"] = np.zeros((natom, 1))
#     j_dict["diagonal"] = np.zeros((natom, 9))
#     j_dict["tensor"] = np.zeros((natom, 9))
#     j_dict["dmi"] = np.zeros((natom, 3))
#     j_dict["symmetric"] = np.zeros((natom, 3))
#     for irow, row in enumerate(exchange):
#         j_mat = row[5:14].reshape(3, 3)
#         j_mat_sym = 0.5 * (j_mat - j_mat.T)
#         j_mat_asym = 0.5 * (j_mat - j_mat.T)
#         dmi = np.array([j_mat_asym[1, 2], j_mat_asym[2, 0], j_mat_asym[0, 1]])
#         jsy = np.array([j_mat_sym[1, 2], j_mat_sym[2, 0], j_mat_sym[0, 1]])
#         j_dict["scalar"][irow] = np.trace(j_mat) / 3.0
#         j_dict["diagonal"][irow] = np.diag(np.diag(j_mat)).flatten()
#         j_dict["dmi"][irow] = dmi
#         j_dict["symmetric"][irow] = jsy
#         j_dict["tensor"][irow] = j_mat.flatten()
# 
#     j_dict["left"] = exchange[:, 0:5]
#     j_dict["right"] = exchange[:, 14:]
#     
#     return j_dict
# 
# def print_projections(j_dict):
# 
#     fmt_l = "%4d %4d   % 4.1f % 4.1f % 4.1f    "
#     fmt_1 = "% 10.6f"
#     fmt_3 = "% 10.6f % 10.6f % 10.6f"
#     fmt_9 = "% 10.6f % 10.6f % 10.6f  % 10.6f % 10.6f % 10.6f  % 10.6f % 10.6f % 10.6f"
#     fmt_r = "     %8.4f"
#     
#     keys = ["scalar", "diagonal", "dmi", "symmetric", "tensor"]
#     fmts = [fmt_1, fmt_9, fmt_3, fmt_3, fmt_9]
#     
#     for idx, key in enumerate(keys):
#         outmat = np.hstack((j_dict["left"], j_dict[key], j_dict["right"]))
#         fname = f"j_{key}.dat"
#         np.savetxt(fname, outmat, fmt=fmt_l + fmts[idx] + fmt_r)
#     print("Projections saved to files.")
#     return


# List of input files
input_directions = ["100", "010", "001"]
input_atoms = ["1", "2", "3", "4"]

input_files = []

exchange_data = {}

# Extract information from each input file
for direction in input_directions:
    exchange_data[direction] = []
    for atoms in input_atoms:
        file_name = f"spin-{direction}/out-{atoms}"
        extracted_data = rspt.RsptExchange(file_name)
        exchange_data[direction].append(extracted_data.outmap)

concatenated_data = {}
for key, value in exchange_data.items():
    concatenated_data[key] = np.concatenate(value)

x_mask = np.array([[0.0, 0.0, 0.0], [0.0, 0.5, 1.0], [0.0, 1.0, 0.5]])
y_mask = np.array([[0.5, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.5]])
z_mask = np.array([[0.5, 1.0, 0.0], [1.0, 0.5, 0.0], [0.0, 0.0, 0.0]])
mask_list = [x_mask, y_mask, z_mask]

masked_data = {key: value.copy() for key, value in concatenated_data.items()}

for idx, (key, value) in enumerate(masked_data.items()):
    for row in value:
        j_mat = row[5:14].reshape(3, 3)
        j_mat = mask_list[idx] * j_mat * 3.0
        row[5:14] = j_mat.flatten()


summed_data = (
    np.array(masked_data["001"] + masked_data["010"] + masked_data["100"]) / 3.0
)

# fmt1 = "%4d %4d   % 4.1f % 4.1f % 4.1f   "
# fmt2 = "% 10.6f % 10.6f % 10.6f  % 10.6f % 10.6f % 10.6f  % 10.6f % 10.6f % 10.6f     %8.4f"
# fmt = fmt1 + fmt2
# np.savetxt("jtot.txt", summed_data, fmt=fmt)

j_dict = rspt.extract_projections(summed_data)
rspt.print_projections(j_dict)
print("Extraction and storage completed successfully!")
