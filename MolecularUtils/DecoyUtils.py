import random
from typing import Dict, List, Optional, Tuple
from difflib import SequenceMatcher

class DecoyUtils:
    """
    用于生成蛋白质或肽段的 decoy 序列的工具类
    
    支持:
    1. reverse 或 shuffle 方法生成 decoy
    2. 可选择保留末端氨基酸不变
    3. 可过滤与 target 相似度高的 decoy
    """
    
    @staticmethod
    def generate_decoy(sequence: str, 
                       modifications: Optional[Dict[int, str]] = None, 
                       method: str = "reverse", 
                       keep_terminals: bool = True,
                       similarity_threshold: float = 0.5,
                       max_attempts: int = 10) -> Tuple[str, Dict[int, str]]:
        """
        生成 decoy 序列
        
        参数:
            sequence: 原始序列
            modifications: 修饰信息，格式为 {位置: 修饰名称}
            method: 生成方法，"reverse" 或 "shuffle"
            keep_terminals: 是否保留末端氨基酸不变
            similarity_threshold: 相似度阈值，高于此值的 decoy 将被过滤
            max_attempts: 最大尝试次数，用于生成低相似度的 decoy
            
        返回:
            decoy 序列和对应的修饰信息
        """
        if not sequence:
            return "", {}
            
        if modifications is None:
            modifications = {}
            
        if method not in ["reverse", "shuffle"]:
            raise ValueError("method 必须为 'reverse' 或 'shuffle'")
            
        # 尝试生成满足相似度要求的 decoy
        for _ in range(max_attempts):
            if method == "reverse":
                decoy_seq, decoy_mods = DecoyUtils._reverse_sequence(sequence, modifications, keep_terminals)
            else:  # shuffle
                decoy_seq, decoy_mods = DecoyUtils._shuffle_sequence(sequence, modifications, keep_terminals)
                
            # 检查相似度
            similarity = DecoyUtils.calculate_similarity(sequence, decoy_seq)
            if similarity <= similarity_threshold:
                return decoy_seq, decoy_mods
                
        # 如果多次尝试后仍无法满足相似度要求，返回空字符串和空字典
        return "", {}
    
    @staticmethod
    def _reverse_sequence(sequence: str, 
                         modifications: Dict[int, str], 
                         keep_terminals: bool) -> Tuple[str, Dict[int, str]]:
        """反转序列并更新修饰位置"""
        if len(sequence) <= 2:
            return sequence, modifications
            
        if keep_terminals:
            # 保留首尾氨基酸，反转中间部分
            middle = sequence[1:-1]
            reversed_middle = middle[::-1]
            decoy_seq = sequence[0] + reversed_middle + sequence[-1]
            
            # 更新修饰位置
            decoy_mods = {}
            for pos, mod in modifications.items():
                if pos == 0 or pos == len(sequence) - 1:
                    # 首尾位置保持不变
                    decoy_mods[pos] = mod
                else:
                    # 中间位置需要重新映射
                    # 原位置在中间部分的相对位置
                    relative_pos = pos - 1
                    # 反转后的相对位置
                    new_relative_pos = len(middle) - 1 - relative_pos
                    # 转换回绝对位置
                    new_pos = new_relative_pos + 1
                    decoy_mods[new_pos] = mod
        else:
            # 完全反转序列
            decoy_seq = sequence[::-1]
            
            # 更新修饰位置
            decoy_mods = {}
            for pos, mod in modifications.items():
                new_pos = len(sequence) - 1 - pos
                decoy_mods[new_pos] = mod
                
        return decoy_seq, decoy_mods
    
    @staticmethod
    def _shuffle_sequence(sequence: str, 
                         modifications: Dict[int, str], 
                         keep_terminals: bool) -> Tuple[str, Dict[int, str]]:
        """打乱序列并更新修饰位置"""
        if len(sequence) <= 2:
            return sequence, modifications
            
        # 创建位置到氨基酸的映射
        pos_to_aa = {i: aa for i, aa in enumerate(sequence)}
        
        if keep_terminals:
            # 保留首尾氨基酸，打乱中间部分
            middle_pos = list(range(1, len(sequence) - 1))
            shuffled_middle_pos = middle_pos.copy()
            random.shuffle(shuffled_middle_pos)
            
            # 创建原位置到新位置的映射
            pos_mapping = {0: 0, len(sequence) - 1: len(sequence) - 1}  # 首尾位置不变
            for old_pos, new_pos in zip(middle_pos, shuffled_middle_pos):
                pos_mapping[old_pos] = new_pos
                
            # 构建 decoy 序列
            decoy_chars = [''] * len(sequence)
            for old_pos, aa in pos_to_aa.items():
                new_pos = pos_mapping[old_pos]
                decoy_chars[new_pos] = aa
            decoy_seq = ''.join(decoy_chars)
            
            # 更新修饰位置
            decoy_mods = {}
            for pos, mod in modifications.items():
                decoy_mods[pos_mapping[pos]] = mod
        else:
            # 完全打乱序列
            positions = list(range(len(sequence)))
            shuffled_positions = positions.copy()
            random.shuffle(shuffled_positions)
            
            # 创建原位置到新位置的映射
            pos_mapping = {old_pos: new_pos for old_pos, new_pos in zip(positions, shuffled_positions)}
            
            # 构建 decoy 序列
            decoy_chars = [''] * len(sequence)
            for old_pos, aa in pos_to_aa.items():
                new_pos = pos_mapping[old_pos]
                decoy_chars[new_pos] = aa
            decoy_seq = ''.join(decoy_chars)
            
            # 更新修饰位置
            decoy_mods = {}
            for pos, mod in modifications.items():
                decoy_mods[pos_mapping[pos]] = mod
                
        return decoy_seq, decoy_mods
    
    @staticmethod
    def calculate_similarity(seq1: str, seq2: str) -> float:
        """计算两个序列的相似度"""
        if not seq1 or not seq2:
            return 0.0
        try:
            matcher = SequenceMatcher(None, seq1, seq2)
            return matcher.ratio()
        except Exception:
            return 0.0
    
    @staticmethod
    def generate_decoy_batch(sequences: List[str], 
                            modifications_list: Optional[List[Dict[int, str]]] = None, 
                            **kwargs) -> List[Tuple[str, Dict[int, str]]]:
        """
        批量生成 decoy 序列
        
        参数:
            sequences: 原始序列列表
            modifications_list: 修饰信息列表
            **kwargs: 传递给 generate_decoy 的其他参数
            
        返回:
            decoy 序列和对应修饰信息的列表
        """
        if modifications_list is None:
            modifications_list = [{}] * len(sequences)
            
        if len(sequences) != len(modifications_list):
            raise ValueError("sequences 和 modifications_list 长度必须相同")
            
        results = []
        for seq, mods in zip(sequences, modifications_list):
            decoy_seq, decoy_mods = DecoyUtils.generate_decoy(seq, mods, **kwargs)
            results.append((decoy_seq, decoy_mods))
            
        return results
    

if __name__ == "__main__":
    # 示例用法
    sequence = "PEPTIDEK"
    modifications = {2: "Phospho", 5: "Oxidation"}
    
    # 生成 reverse decoy
    decoy_seq, decoy_mods = DecoyUtils.generate_decoy(
        sequence,
        modifications,
        method="reverse", 
        keep_terminals=True
    )
    print(f"原始序列: {sequence}")
    print(f"修饰: {modifications}")
    print(f"Reverse Decoy: {decoy_seq}")
    print(f"Decoy 修饰: {decoy_mods}")
    
    # 生成 shuffle decoy
    decoy_seq, decoy_mods = DecoyUtils.generate_decoy(
        sequence,
        modifications,
        method="shuffle",
        keep_terminals=True
    )
    print(f"\nShuffle Decoy: {decoy_seq}")
    print(f"Decoy 修饰: {decoy_mods}")

    # 解析带修饰标记的序列
    modified_sequence = "PEP(Phospho)TI(Oxidation)DEK"
    clean_seq, mods = DecoyUtils.parse_modified_sequence(modified_sequence)
    print(f"\n带修饰序列: {modified_sequence}")
    print(f"解析后序列: {clean_seq}")
    print(f"解析后修饰: {mods}")

    # 格式化带修饰的序列
    formatted = DecoyUtils.format_modified_sequence(clean_seq, mods)
    print(f"格式化后: {formatted}")

    