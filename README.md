# OxStressLAB-RNAseq: Пайплайн для оценки влияния окислительного стресса на дифференциальную экспрессию генов у молочнокислых бактерий.

![Bioconda](https://img.shields.io/badge/Built_with-R%20%7C%20Bioconductor-blue)
![License](https://img.shields.io/badge/License-MIT-green)

## Описание
Пайплайн для анализа дифференциальной экспрессии генов у молочнокислых бактерий (LAB) при окислительном стрессе. Включает:
- Контроль качества данных
- Предобработку и нормализацию
- Анализ дифференциальной экспрессии (DESeq2)
- Визуализацию результатов
- Функциональный анализ (GO/KEGG)

## Системные требования
| Компонент | Минимальные требования |
|-----------|------------------------|
| ОС | Windows 10+/macOS 10.15+/Linux Ubuntu 20.04+ |
| R | версия 4.2.0 или новее |
| RStudio | 2023.03+ (рекомендуется) |
| Память | 8 ГБ ОЗУ (16+ ГБ для больших датасетов) |
| Место на диске | 5 ГБ свободного пространства |

## Установка

### 1. Установка R и RStudio
```bash
# Для Linux (Ubuntu/Debian):
sudo apt-get install r-base

## Установка пакетов
В R-консоли выполните:
### Установка менеджера пакетов Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

### Основные пакеты
BiocManager::install(c(
    "DESeq2",       # Анализ дифференциальной экспрессии
    "apeglm",       # Уточнение оценок log2FC
    "pheatmap",     # Визуализация тепловых карт
    "RColorBrewer", # Цветовые палитры
    "clusterProfiler" # Функциональный анализ
))

### Проверка установки
library(DESeq2)

## Запуск пайплайна
1. Подготовка данных
Структура файлов:

project/
├── data/
│   ├── counts.csv     # Таблица счетов (гены × образцы)
│   └── metadata.csv   # Метаданные экспериментов
└── scripts/
    └── rna_seq_pipeline.R
Пример counts.csv:

csv
gene,control_1,control_2,stress_1,stress_2
gene1,150,200,50,30
gene2,3000,2800,1000,900
2. Выполнение анализа
r
### Установка рабочей директории
setwd("/path/to/project")

### Запуск пайплайна
source("scripts/rna_seq_pipeline.R")

## Частые ошибки и решения
### 1. Ошибка установки пакетов
r
Error: package 'DESeq2' not found
Решение:

r
### Обновите Bioconductor
BiocManager::install(version = "3.18")

### Принудительная установка
BiocManager::install("DESeq2", force = TRUE)

### 2. Отсутствие apeglm
r
Error in lfcShrink(..., type="apeglm"): 
  type='apeglm' requires the Bioconductor package 'apeglm'
Решение:

r
BiocManager::install("apeglm")

### 3. Проблема с vst()
r
Warning: less than 'nsub' rows in count matrix
Альтернатива:

r
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

### 4. Проблемы с памятью
Симптомы:

RStudio зависает

Сообщения "cannot allocate vector of size..."

Решение:

r
### Увеличьте лимит памяти
options(future.globals.maxSize = 8000 * 1024^2) # 8GB

## Поддержка
При возникновении проблем:

Проверьте лог-файл logs/analysis.log

Соберите информацию о системе:

r
sessionInfo()

Создайте Issue в репозитории с:
- Текст ошибки
- Пример входных данных
- Вывод sessionInfo()
