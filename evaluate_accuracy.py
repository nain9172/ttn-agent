#!/usr/bin/env python3
"""
Evaluate TTN Agent Accuracy using ClinVar Dataset
評估 TTN Agent 使用 ClinVar 資料集的精準度
"""

import argparse
import csv
import sys
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple
import json

from main import setup_logging
from utils.variant_parser import parse_variant
from utils.evo2_predictor import Evo2Predictor
from utils.clinvar_parser import ClinVarParser
from config import OUTPUT_DIR


class AccuracyEvaluator:
    """評估 TTN Agent 預測精準度"""
    
    def __init__(self, skip_evo2: bool = False):
        self.logger = logging.getLogger(__name__)
        self.skip_evo2 = skip_evo2
        self.evo2_predictor = None if skip_evo2 else Evo2Predictor()
        self.clinvar_parser = ClinVarParser()
        
        # 統計數據
        self.stats = {
            'total': 0,
            'processed': 0,
            'failed': 0,
            'evo2_predictions': {
                'pathogenic': 0,
                'benign': 0,
                'uncertain': 0
            },
            'clinvar_labels': {
                'pathogenic': 0,
                'likely_pathogenic': 0,
                'benign': 0,
                'likely_benign': 0,
                'uncertain': 0,
                'conflicting': 0
            },
            'confusion_matrix': {
                'true_positive': 0,   # 正確預測為致病
                'true_negative': 0,   # 正確預測為良性
                'false_positive': 0,  # 錯誤預測為致病（實際良性）
                'false_negative': 0   # 錯誤預測為良性（實際致病）
            }
        }
        
        # 詳細結果
        self.detailed_results = []
    
    def read_clinvar_csv(self, csv_path: Path) -> List[Dict]:
        """讀取 ClinVar CSV 文件"""
        variants = []
        
        with open(csv_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                variants.append({
                    'classification': row['Germline classification'],
                    'chr': row['chr'],
                    'position': row['position'],
                    'ref': row['ref'],
                    'alt': row['alt']
                })
        
        self.logger.info(f"從 {csv_path} 讀取了 {len(variants)} 個變異")
        return variants
    
    def classify_clinvar_label(self, classification: str) -> str:
        """
        將 ClinVar 分類標準化為 pathogenic/benign/uncertain
        
        分類規則：
        - Pathogenic, Likely pathogenic, Pathogenic/Likely pathogenic → pathogenic
        - Benign, Likely benign, Benign/Likely benign → benign
        - Uncertain significance, Conflicting interpretations → uncertain
        """
        classification_lower = classification.lower()
        
        if 'pathogenic' in classification_lower and 'benign' not in classification_lower:
            return 'pathogenic'
        elif 'benign' in classification_lower:
            return 'benign'
        else:
            return 'uncertain'
    
    def evaluate_variant(self, variant_data: Dict) -> Dict:
        """評估單個變異"""
        try:
            # 構建變異字符串
            variant_str = f"{variant_data['chr']}-{variant_data['position']}-{variant_data['ref']}-{variant_data['alt']}"
            
            # 解析變異
            variant_info = parse_variant(variant_str)
            
            # ClinVar 分類（ground truth）
            clinvar_classification = variant_data['classification']
            clinvar_label = self.classify_clinvar_label(clinvar_classification)
            
            # Evo2 預測
            evo2_prediction = None
            evo2_delta = None
            evo2_confidence = None
            
            if not self.skip_evo2:
                evo2_result = self.evo2_predictor.predict(variant_info)
                if evo2_result and evo2_result.get('success'):
                    evo2_prediction = evo2_result.get('prediction')
                    evo2_delta = evo2_result.get('delta_score')
                    evo2_confidence = evo2_result.get('confidence')
            
            # ClinVar API 查詢（可選）
            clinvar_info = self.clinvar_parser.parse_variant(variant_info)
            clinvar_api_significance = None
            if clinvar_info:
                clinvar_api_significance = clinvar_info.get('clinical_significance', 'Not found')
            
            result = {
                'variant_id': variant_info['variant_id'],
                'clinvar_classification': clinvar_classification,
                'clinvar_label': clinvar_label,
                'evo2_prediction': evo2_prediction,
                'evo2_delta': evo2_delta,
                'evo2_confidence': evo2_confidence,
                'clinvar_api_significance': clinvar_api_significance,
                'success': True
            }
            
            # 更新統計
            self.update_statistics(result)
            
            return result
            
        except Exception as e:
            self.logger.error(f"評估變異時出錯 {variant_data}: {e}")
            return {
                'variant_id': f"{variant_data['chr']}-{variant_data['position']}-{variant_data['ref']}-{variant_data['alt']}",
                'clinvar_classification': variant_data['classification'],
                'clinvar_label': self.classify_clinvar_label(variant_data['classification']),
                'evo2_prediction': None,
                'evo2_delta': None,
                'evo2_confidence': None,
                'clinvar_api_significance': None,
                'success': False,
                'error': str(e)
            }
    
    def update_statistics(self, result: Dict):
        """更新統計數據"""
        self.stats['processed'] += 1
        
        # ClinVar 標籤統計
        clinvar_label = result['clinvar_label']
        if clinvar_label == 'pathogenic':
            self.stats['clinvar_labels']['pathogenic'] += 1
        elif clinvar_label == 'benign':
            self.stats['clinvar_labels']['benign'] += 1
        else:
            self.stats['clinvar_labels']['uncertain'] += 1
        
        # Evo2 預測統計
        evo2_pred = result['evo2_prediction']
        if evo2_pred:
            if evo2_pred == 'pathogenic':
                self.stats['evo2_predictions']['pathogenic'] += 1
            elif evo2_pred == 'benign':
                self.stats['evo2_predictions']['benign'] += 1
            else:
                self.stats['evo2_predictions']['uncertain'] += 1
            
            # 更新混淆矩陣（只對確定的標籤進行評估）
            if clinvar_label in ['pathogenic', 'benign'] and evo2_pred in ['pathogenic', 'benign']:
                if clinvar_label == 'pathogenic' and evo2_pred == 'pathogenic':
                    self.stats['confusion_matrix']['true_positive'] += 1
                elif clinvar_label == 'benign' and evo2_pred == 'benign':
                    self.stats['confusion_matrix']['true_negative'] += 1
                elif clinvar_label == 'benign' and evo2_pred == 'pathogenic':
                    self.stats['confusion_matrix']['false_positive'] += 1
                elif clinvar_label == 'pathogenic' and evo2_pred == 'benign':
                    self.stats['confusion_matrix']['false_negative'] += 1
    
    def calculate_metrics(self) -> Dict:
        """計算評估指標"""
        cm = self.stats['confusion_matrix']
        tp = cm['true_positive']
        tn = cm['true_negative']
        fp = cm['false_positive']
        fn = cm['false_negative']
        
        # 準確率 (Accuracy)
        total_evaluated = tp + tn + fp + fn
        accuracy = (tp + tn) / total_evaluated if total_evaluated > 0 else 0
        
        # 精確率 (Precision) - 預測為致病的準確性
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        
        # 召回率 (Recall/Sensitivity) - 找出致病變異的能力
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        
        # 特異性 (Specificity) - 找出良性變異的能力
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        
        # F1 分數
        f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        
        return {
            'total_variants': self.stats['total'],
            'processed': self.stats['processed'],
            'failed': self.stats['failed'],
            'evaluated': total_evaluated,
            'accuracy': accuracy,
            'precision': precision,
            'recall': recall,
            'specificity': specificity,
            'f1_score': f1_score,
            'confusion_matrix': cm
        }
    
    def evaluate_dataset(self, csv_path: Path) -> Dict:
        """評估整個數據集"""
        self.logger.info(f"開始評估數據集: {csv_path}")
        
        # 讀取變異
        variants = self.read_clinvar_csv(csv_path)
        self.stats['total'] = len(variants)
        
        # print(f"\n{'='*70}")
        # print(f"開始處理 {len(variants)} 個變異")
        # print(f"{'='*70}\n")
        
        # 評估每個變異
        for idx, variant_data in enumerate(variants, 1):
            variant_str = f"{variant_data['chr']}-{variant_data['position']}-{variant_data['ref']}-{variant_data['alt']}"
            
            # 即時顯示當前處理的變異
            # print(f"\r[{idx}/{len(variants)}] 處理中: {variant_str}", end='', flush=True)
            
            # self.logger.info(f"\n{'='*70}")
            # self.logger.info(f"評估變異 {idx}/{len(variants)}: {variant_str}")
            # self.logger.info(f"{'='*70}")
            
            result = self.evaluate_variant(variant_data)
            self.detailed_results.append(result)
            
            if not result['success']:
                self.stats['failed'] += 1
            
            # 每完成一個變異就顯示詳細進度
            metrics = self.calculate_metrics()
            clinvar_label = result['clinvar_label']
            evo2_pred = result.get('evo2_prediction', 'N/A')
            
            # 清除之前的行並顯示完整信息
            # print(f"\r[{idx}/{len(variants)}] {variant_str} | ClinVar: {clinvar_label} | Evo2: {evo2_pred} | 準確率: {metrics['accuracy']:.1%} ({metrics['evaluated']} 個已評估)")
            
            # 每10個顯示詳細統計
            if idx % 10 == 0:
                # print(f"  ├─ 致病: {self.stats['clinvar_labels']['pathogenic']}, 良性: {self.stats['clinvar_labels']['benign']}")
                if metrics['evaluated'] > 0:
                    cm = metrics['confusion_matrix']
        #             print(f"  └─ TP:{cm['true_positive']} TN:{cm['true_negative']} FP:{cm['false_positive']} FN:{cm['false_negative']}")
        
        # print(f"\n{'='*70}")
        # print("評估完成！")
        # print(f"{'='*70}\n")
        
        # 計算最終指標
        final_metrics = self.calculate_metrics()
        
        return final_metrics
    
    def save_results(self, metrics: Dict, output_dir: Path):
        """保存評估結果"""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 保存詳細結果 CSV
        csv_path = output_dir / 'detailed_results.csv'
        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            fieldnames = [
                'variant_id', 'clinvar_classification', 'clinvar_label',
                'evo2_prediction', 'evo2_delta', 'evo2_confidence',
                'clinvar_api_significance', 'success'
            ]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for result in self.detailed_results:
                # 只寫入需要的欄位
                row = {k: result.get(k, '') for k in fieldnames}
                writer.writerow(row)
        
        self.logger.info(f"詳細結果已保存至: {csv_path}")
        
        # 保存統計摘要 JSON
        summary_path = output_dir / 'summary.json'
        summary = {
            'statistics': self.stats,
            'metrics': metrics,
            'evaluation_date': datetime.now().isoformat()
        }
        
        with open(summary_path, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        
        self.logger.info(f"統計摘要已保存至: {summary_path}")
        
        # 保存可讀性報告
        report_path = output_dir / 'evaluation_report.txt'
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("="*70 + "\n")
            f.write("TTN Agent 精準度評估報告\n")
            f.write("="*70 + "\n\n")
            
            f.write("資料集統計\n")
            f.write("-"*70 + "\n")
            f.write(f"總變異數: {metrics['total_variants']}\n")
            f.write(f"成功處理: {metrics['processed']}\n")
            f.write(f"失敗: {metrics['failed']}\n")
            f.write(f"參與評估: {metrics['evaluated']}\n\n")
            
            f.write("ClinVar 標籤分布\n")
            f.write("-"*70 + "\n")
            for label, count in self.stats['clinvar_labels'].items():
                f.write(f"{label}: {count}\n")
            f.write("\n")
            
            f.write("Evo2 預測分布\n")
            f.write("-"*70 + "\n")
            for pred, count in self.stats['evo2_predictions'].items():
                f.write(f"{pred}: {count}\n")
            f.write("\n")
            
            f.write("混淆矩陣\n")
            f.write("-"*70 + "\n")
            cm = metrics['confusion_matrix']
            f.write(f"True Positive (正確預測致病):  {cm['true_positive']}\n")
            f.write(f"True Negative (正確預測良性):  {cm['true_negative']}\n")
            f.write(f"False Positive (誤報致病):     {cm['false_positive']}\n")
            f.write(f"False Negative (漏報致病):     {cm['false_negative']}\n\n")
            
            f.write("評估指標\n")
            f.write("-"*70 + "\n")
            f.write(f"準確率 (Accuracy):    {metrics['accuracy']:.2%}\n")
            f.write(f"精確率 (Precision):   {metrics['precision']:.2%}\n")
            f.write(f"召回率 (Recall):      {metrics['recall']:.2%}\n")
            f.write(f"特異性 (Specificity): {metrics['specificity']:.2%}\n")
            f.write(f"F1 分數:              {metrics['f1_score']:.4f}\n\n")
            
            f.write("="*70 + "\n")
        
        self.logger.info(f"評估報告已保存至: {report_path}")


def main():
    """主函數"""
    parser = argparse.ArgumentParser(
        description='評估 TTN Agent 使用 ClinVar 資料集的精準度'
    )
    parser.add_argument(
        'csv_file',
        type=str,
        default='/home/ryan910702/ttn-agent/data/clinvar_result_with_ref_alt.csv',
        nargs='?',
        help='ClinVar CSV 檔案路徑'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='輸出目錄'
    )
    parser.add_argument(
        '--skip-evo2',
        action='store_true',
        help='跳過 Evo2 預測（僅測試 ClinVar API）'
    )
    parser.add_argument(
        '--limit',
        type=int,
        default=None,
        help='限制評估的變異數量（用於快速測試）'
    )
    
    args = parser.parse_args()
    
    # 設置日誌
    logger = setup_logging()
    
    # 檢查輸入文件
    csv_path = Path(args.csv_file)
    if not csv_path.exists():
        logger.error(f"找不到檔案: {csv_path}")
        return 1
    
    # 設置輸出目錄
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = OUTPUT_DIR / f"evaluation_{timestamp}"
    
    # logger.info(f"輸出目錄: {output_dir}")
    
    # 創建評估器
    evaluator = AccuracyEvaluator(skip_evo2=args.skip_evo2)
    
    # 如果有限制，只讀取部分變異
    if args.limit:
        # logger.info(f"限制評估數量: {args.limit}")
        # 修改讀取邏輯以支持限制
        original_read = evaluator.read_clinvar_csv
        def limited_read(path):
            variants = original_read(path)
            return variants[:args.limit]
        evaluator.read_clinvar_csv = limited_read
    
    # 評估數據集
    # logger.info("\n" + "="*70)
    # logger.info("開始評估 TTN Agent 精準度")
    # logger.info("="*70 + "\n")
    
    metrics = evaluator.evaluate_dataset(csv_path)
    
    # 保存結果
    evaluator.save_results(metrics, output_dir)
    
    # 顯示結果
    # print("\n" + "="*70)
    # print("評估完成！")
    # print("="*70)
    # print(f"\n總變異數:     {metrics['total_variants']}")
    # print(f"成功處理:     {metrics['processed']}")
    # print(f"失敗:         {metrics['failed']}")
    # print(f"參與評估:     {metrics['evaluated']}")
    # print(f"\n準確率:       {metrics['accuracy']:.2%}")
    # print(f"精確率:       {metrics['precision']:.2%}")
    # print(f"召回率:       {metrics['recall']:.2%}")
    # print(f"特異性:       {metrics['specificity']:.2%}")
    # print(f"F1 分數:      {metrics['f1_score']:.4f}")
    # print(f"\n輸出目錄:     {output_dir}")
    # print("="*70 + "\n")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

