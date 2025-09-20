class FDRUtils():
    def __init__(self):
        pass

    def calculate_fdr_list(self, score, label):
        """
        Calculate the False Discovery Rate (FDR) based on the given score and label.

        Parameters:
        score (float): The score to evaluate.
        label (int): The label indicating whether the score is positive (1) or negative (0).

        Returns:
        list: A list containing the FDR values.
        """
        # Sort scores in descending order and track original indices
        sorted_indices = sorted(range(len(score)), key=lambda i: score[i], reverse=True)

        # Create sorted lists of scores and labels
        sorted_scores = [score[i] for i in sorted_indices]
        sorted_labels = [label[i] for i in sorted_indices]

        # Calculate FDR values
        fdr_list = []
        target_count = 0
        decoy_count = 0

        for i in range(len(sorted_scores)):
            if sorted_labels[i] == 1:
                target_count += 1
            else:
                decoy_count += 1

            if target_count == 0:
                fdr_list.append(0)
            else:
                fdr_list.append(float(decoy_count) / target_count)

        # Calculate q-values using monotonic FDR
        min_fdr = fdr_list[-1]
        fdr_list_mono = []

        for i in range(len(fdr_list) - 1, -1, -1):
            if min_fdr > fdr_list[i]:
                min_fdr = fdr_list[i]
            fdr_list_mono.append(min_fdr)

        fdr_list_mono.reverse()

        # Map q-values back to original order
        result_list = [0] * len(score)
        for i, orig_idx in enumerate(sorted_indices):
            result_list[orig_idx] = fdr_list_mono[i]

        return result_list

    def calculate_fdr(self, score, label, target_fdr=0.01, top_n=20):
        """
        Calculate FDR metrics including count of targets below threshold and minimum score threshold.
        Before FDR calculation, removes top 20 decoy (label=0) data points.

        Parameters:
        score (list): List of scores
        label (list): List of labels (1 for target, 0 for decoy)
        target_fdr (float): Target FDR threshold (default: 0.01)

        Returns:
        tuple: (count_below_threshold, min_score_threshold)
            - count_below_threshold: Number of targets with FDR below threshold
            - min_score_threshold: Minimum score that satisfies target FDR
        """
        # Create list of (score, label) pairs
        score_label_pairs = list(zip(score, label))
        
        # Sort by score in descending order
        score_label_pairs.sort(key=lambda x: x[0], reverse=True)
        
        # Remove top top_n decoy data points
        decoy_count = 0
        filtered_pairs = []
        for s, l in score_label_pairs:
            if l == 0 and decoy_count < top_n:
                decoy_count += 1
                continue
            filtered_pairs.append((s, l))
            
        # Unzip filtered pairs back into separate lists
        filtered_score, filtered_label = zip(*filtered_pairs)
        
        # Calculate FDR on filtered data
        fdr_list = self.calculate_fdr_list(filtered_score, filtered_label)
        count_below_threshold = sum(1 for index, fdr in enumerate(fdr_list) if fdr < target_fdr and filtered_label[index] == 1)
        
        min_score = float('inf')
        for i, fdr in enumerate(fdr_list):
            if fdr <= target_fdr:
                min_score = min(min_score, filtered_score[i])
                
        return count_below_threshold, min_score
    