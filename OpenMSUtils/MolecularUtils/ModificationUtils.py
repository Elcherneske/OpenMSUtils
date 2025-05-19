from .accurate_molmass import EnhancedFormula
import os
import pandas as pd
from typing import Tuple, Dict, List, Optional
import re

class Modification():
    def __init__(self, name: str, formula: str):
        self.name = name
        self.formula = formula

    @property
    def mass(self):
        # 处理特殊的加合物表示法，如 [M+H+], [M+NH3] 等
        adduct_pattern = r'\[M([\+\-].+)\]'
        match = re.match(adduct_pattern, self.formula)
        if match:
            adduct = match.group(1)
            if adduct.startswith('+'):
                return EnhancedFormula(adduct[1:]).isotope.mass
            elif adduct.startswith('-'):
                return -EnhancedFormula(adduct[1:]).isotope.mass
            else:
                raise ValueError(f"Invalid adduct format: {self.formula}")
        else:
            raise ValueError(f"Invalid adduct format: {self.formula}")

class ModificationUtils():
    @staticmethod
    def parse_modification_file(file_path: str) -> pd.DataFrame:
        """
        加载修饰库文件
        
        参数:
            file_path: 修饰库文件路径
        """
        if os.path.exists(file_path):
            return pd.read_csv(file_path, sep='\t')
        else:
            raise FileNotFoundError(f"File {file_path} not found")
    
    @staticmethod
    def find_modifications_by_mass(modifications: pd.DataFrame, mass: float, tolerance: float = 0.0001) -> List[Modification]:
        """
        根据质量搜索修饰
        
        参数:
            mass: 修饰质量
            tolerance: 质量容差
            
        返回:
            符合条件的修饰列表
        """
        if modifications is None:
            return []
            
        results = []
        for modification in modifications:
            if abs(modification['Isotopic Mass'] - mass) <= tolerance:
                results.append(Modification(modification['Mod Name'], modification['Formula']))
        return results
    
    @staticmethod
    def parse_modified_sequence(modified_sequence: str) -> Tuple[str, Dict[int, str]]:
        """
        解析带修饰标记的序列，支持多种修饰格式：
        - PEP(Phospho)TIDE
        - PEP(Phospho (P))TIDE
        - PEP[Phospho (P)]TIDE(Phospho)
        - PEP[Phospho [P]]TIDE[Phospho (P)]
        
        参数:
            modified_sequence: 带修饰标记的序列
            
        返回:
            原始序列和修饰信息
        """
        clean_sequence = ""
        modifications = {}
        i = 0
        pos = 0
        
        while i < len(modified_sequence):
            if i < len(modified_sequence) - 1 and modified_sequence[i+1] in ['(', '[']:
                # 找到一个氨基酸后面跟着修饰
                aa = modified_sequence[i]
                clean_sequence += aa
                pos = len(clean_sequence) - 1
                
                # 确定修饰的开始和结束符号
                start_symbol = modified_sequence[i+1]
                end_symbol = ')' if start_symbol == '(' else ']'
                
                # 找到完整的修饰（处理嵌套括号）
                start = i + 2  # 跳过开始符号
                paren_count = 1
                j = start
                
                while j < len(modified_sequence) and paren_count > 0:
                    if modified_sequence[j] in ['(', '[']:
                        paren_count += 1
                    elif modified_sequence[j] in [')', ']']:
                        paren_count -= 1
                    j += 1
                
                if paren_count == 0:
                    mod = modified_sequence[start:j-1]  # 去除最外层括号
                    modifications[pos] = mod
                    i = j  # 更新索引到修饰后的位置
                    continue
            
            # 普通字符
            clean_sequence += modified_sequence[i]
            i += 1
        
        return clean_sequence, modifications
    
    @staticmethod
    def format_modified_sequence(sequence: str, modifications: Dict[int, str]) -> str:
        """
        将序列和修饰信息格式化为带修饰标记的序列
        
        参数:
            sequence: 原始序列
            modifications: 修饰信息
            
        返回:
            带修饰标记的序列
        """
        result = ""
        for i, aa in enumerate(sequence):
            result += aa
            if i in modifications:
                result += f"({modifications[i]})"
                
        return result


