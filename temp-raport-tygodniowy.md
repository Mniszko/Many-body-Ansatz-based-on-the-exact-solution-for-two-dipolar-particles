1. Zmieniłem program tak, żeby możliwe było generowanie danych na temat dyspersji
    - występują problemy z całkowaniem - przy niektórych wartościach parametrów (sprawdzałem długość) dyspersja jest zespolona, co jest niemożliwe
    - całkowite całki po $r$ i $r^2$ są znacznie większe w tych przypadkach, chociaż co jest oczywiste być nie powinny
    - wraz z odległością jeśli całka jest liczona (wydaje się) poprawnie dyspersja dąży do zera, co też jest bardzo niepokojące
2. zapisuję dane w postaci plików csv - nie macierze, ale chociażby wektory własne systemu już tak
3. napisałem kod do liczenia oddziaływania dipolowego
    - można go przyspieszyć uwzględniając specjalne wartości liczb kwantowych, np. zmniejszyć rząd o przynajmniej 1 kiedy n1=n3 i tak dalej
    - narazie całkowanie przez dowolny przedział wywala się kompletnie, sprubuję wykonać podobną całkę w kontrolowanym środowisku, być może z podstawieniem które sam sobie wymyślę.



# CODE REVIEWER

no więc jest problem natury posikanej
norm1 i norm2 to współczynniki normalizacyjne pojedyńczych funkcji i są przetrzymywane w liście w parach. teraz. współczynnik normalizacyjny to 1/sqrt(norm)

dodatkowe notatki na stronach 119,120 w zeszycie

Prawdopodobnie mój błąd polegał na tym, że stosowałem normalizację odpowiednią dla dwóch cząstek, w przypadku dyspersji nie może być mowy o normalizacji w wyniku której otrzymujemy dwukrotność prawdopodobieństwa - wtedy wzór na nieoznaczoność obserwabli nie jest poprawny. Trzeba wprowadzić poprawkę normalizacyjną 1/2 - w tym celu określam normsym równą 1/2 albo 1/4 - to jest bardzo ważne, bo przy błędnej normalizacji r2 jest dwukrotnie większe niż powinno a r^2 aż 4 krotnie! Z tego powodu relacja trójkąta dla obserwabli nie jest spełniona

jak już wcześniej wspominałem powinienem dokładnie przebadać jaki iloczyn faktycznie liczyć, być może nie wychwyciłem tego problemu wcześniej, ale całkowanie po całej przestrzeni

# CO MOGĘ ZROBIĆ W SPRAWIE DIPOLI

mógłbym sprubować zastosować pow(x,-5) zamiast 1/x^5, może to coś zmieni, może nie
mógłbym też spróbować zamienić zmienne


wartość dyspersji dąży do 0.3014 - dyspersji u1=0,u2=0,m1=0,m2=0 (tak się spodziewam, faktycznie przy dużych odległościach występują też spore prawdopodobieństwa obsadzenia stanów m (chyba))
