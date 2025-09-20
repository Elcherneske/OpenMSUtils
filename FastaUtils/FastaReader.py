class FastaReader:
    def __init__(self):
        pass

    def read(self, filename):
        with open(filename, 'r') as file:
            lines = file.readlines()
            sequences = self._parse_fasta(lines)
        return sequences

    def _parse_fasta(self, lines):
        """
        line: [str1, str2, ...]
        """
        sequences = {}
        current_header = None
        current_sequence = []
        for line in lines:
            line = line.strip()
            if line.startswith('>'):  # 以 '>' 开头的是 FASTA 头部行
                if current_header:  # 如果已经有一个序列，保存之前的序列
                    sequences[current_header] = ''.join(current_sequence)
                current_header = line[1:]  # 去掉 '>'
                current_sequence = []  # 重置序列
            else:
                current_sequence.append(line)  # 收集序列的每一行
        # 保存最后一个序列
        if current_header:
            sequences[current_header] = ''.join(current_sequence)

        return sequences