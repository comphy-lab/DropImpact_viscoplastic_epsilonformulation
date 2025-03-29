import numpy as np

# 生成 Js 和 Wes
Js = np.append(0, np.logspace(-2, 1, 31))  # 32 个元素
Wes = np.logspace(0, 3, 32)               # 32 个元素

# 分别取前 16 个、后 16 个，并使用 np.round(…, 4) 进行四舍五入保留 4 位小数
# 注意：下面依然会在输出时用 f"{value:.4f}" 强制保留四位小数，以便包括 trailing zeros
Js_first_16  = Js[:16]
Js_second_16 = Js[16:]
Wes_first_16  = Wes[:16]
Wes_second_16 = Wes[16:]

def format_array(array):
    """
    接收一个浮点数数组，返回每个元素都严格以4位小数形式输出的字符串列表
    """
    return [f"{val:.4f}" for val in array]

# 将结果输出到文本文件 output.txt
with open("output.txt", "w") as f:
    # 写入 Js_first_16
    f.write("Js_first_16: ")
    f.write(" ".join(format_array(Js_first_16)))
    f.write("\n")

    # 写入 Js_second_16
    f.write("Js_second_16: ")
    f.write(" ".join(format_array(Js_second_16)))
    f.write("\n")

    # 写入 Wes_first_16
    f.write("Wes_first_16: ")
    f.write(" ".join(format_array(Wes_first_16)))
    f.write("\n")

    # 写入 Wes_second_16
    f.write("Wes_second_16: ")
    f.write(" ".join(format_array(Wes_second_16)))
    f.write("\n")
