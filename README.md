# TRS-omix: search engine

This code accompanies the paper:
Sebastian Sakowski, Marta Majchrzak, Jacek Waldmajer, Pawel Parniewski: TRS-omix: a new search engine for trinucleotide flanked sequences. 2021.

The content of the repository comes from it's predecessor, i.e.:
https://github.com/TRS-omix/software

# ToDo's:
"Proste tematy"
0. Korekta sekwencji flankującej
*CGACGACGACG* + analogicznie po prawej strony

1. Przy uruchomieniu programu plik wsadowy z listą genomów, nazwy plików genomowych dowolne

2. W kolumnie O ("FRG NO") - wpis generowany ręcznie, do automatyzacji.
ciąg znaków z nagłowka pliku fasta (>..." ")  prefix-indeks

3. Podobieństwo sekwencji w ramach iteriors.txt ("wnętrza")

# Opisane propozycje znajdują się w pliku Proposed code.ipynb

1. Program pyta użytkownika o lokalizację plików .fasta i odczytuje wszystkie znajdujące się w danym folderze niezależnie od nazwy - ta metoda wydaję się być odrobinę szybsza
2. Program pyta użytkownika o minimalną oraz maksymalną długość sekwencji jak i również o tryb w którym ma być uruchomiony TRS_omix
3. Połączone pliki interiors.txt są przechowywane jako DataFrame oraz zapisane do pliku .csv
4. Program pyta użytkownika o długość sekwencji którą chcę otrzymać zarówno z początku jak i końca sekwencji
5. Przy użyciu Entrez z Biopython kolumna GENOME służy do przeszukania bazy nucleotide NCBI i znalezienia nazw organizmów
6. Nazwa organizmu,kierunek z którego pochodzi sekwencja(L/R),L-No/R-No oraz indeks sekwencji są wykorzystywane do stworzenia unikalnych identyfikatorów sekwencji
7. Uzyskane w ten sposób pliki .fasta clustrujemy przy pomocy cd-hit-est(https://github.com/weizhongli/cdhit) tak aby odnaleźć sekwencje które pokrywają się w 100% usuwamy je z pierwotnego zestawu sekwencji i przenosimy do innego pliku
8. Blastujemy pozostałe  

# Inne zmiany 
1. environment.yml jest teraz skonfigurowany prawidłowo i umożliwia łatwe odtworzenie środowiska przy pomocy condy
2. dodano brakujące __init__.py
