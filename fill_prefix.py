import pandas as pd
import re

# 读取Excel文件
file_path = './上线货号待剔除.xlsx'
xls = pd.ExcelFile(file_path)
output_xls = "./上线货号剔除后.xlsx"

# 定义处理逻辑
def process_sheet(df: pd.DataFrame, header: list) -> pd.DataFrame:
    # 根据第一张表的表头处理
    for col in df.columns:
        if header[col] == 'i5-1.5':
            df[col] = df[col].fillna('')  # 用空字符串替代NaN
            df[col] = df[col].apply(lambda x: 'AC' + x if x else x)  # 只有非空才添加前缀
        elif header[col] == 'i7-1.0':
            df[col] = df[col].fillna('')  # 用空字符串替代NaN
            df[col] = df[col].apply(lambda x: 'AT' + x if x else x)  # 只有非空才添加前缀
    return df

# 第一张表头
header = pd.read_excel(xls, sheet_name=xls.sheet_names[0], header=None).iloc[0].to_list()
output_table = pd.ExcelWriter(output_xls)
# 处理以"N"或"T"开头并且后续只包含数字和"-"的子表
for sheet_name in xls.sheet_names:
    df = pd.read_excel(xls, sheet_name=sheet_name, header=None)
    if re.match(r'[NT][\d-]+$', sheet_name):
        processed_df = process_sheet(df, header)
        # 输出处理后的表格到Excel文件中
        processed_df.to_excel(output_table, index=False, sheet_name=sheet_name, header=None)
    else:
        df.to_excel(output_table, index=False, sheet_name=sheet_name, header=None)
output_table.close()
print("子表处理完成并输出为Excel文件。")
