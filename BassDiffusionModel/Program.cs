using Microsoft.Research.Oslo;
using System;
using System.Linq;
using System.IO;

namespace BassModel
{
    class Program
    {
        //DO ZMIANY W RAZIE POTRZEBY
        private const string _path = @"bass.dat";

        string[] stringsToRemove = new string[] { "[", "]" };
        private static int[] finalAdopters, initialAdopters;
        private static double[] p, q, u, finalTime, initialTime;
        private static int _setNumber;

        //p - efekt reklam
        //q - wg wikipedii - efekt interakcji między adaptatorami "mouth to mouth"
        //u - wsp adaptatorów, którzy przestają nimi być - wpływ na ujęcie od ostatecznej ilości 
        //initialAdopters - liczba aktualnych kupców
        public static void Start()
        {
            finalAdopters = new int[4];
            initialAdopters = new int[4];
            p = new double[4];
            q = new double[4];
            u = new double[4];
            initialTime = new double[4];
            finalTime = new double[4];

            p[0] = 0.03;
            q[0] = 0.38;
            u[0] = 0.0;
            finalAdopters[0] = 20000;
            initialAdopters[0] = 1;
            initialTime[0] = 0.0;
            finalTime[0] = 100;

            p[1] = 0.00;
            q[1] = 0.5;
            u[1] = 0.0;
            finalAdopters[1] = 20000;
            initialAdopters[1] = 1;
            initialTime[1] = 0.0;
            finalTime[1] = 100;

            p[2] = 0.05;
            q[2] = 0.5;
            u[2] = 0.1;
            finalAdopters[2] = 20000;
            initialAdopters[2] = 10000;
            initialTime[2] = 0.0;
            finalTime[2] = 100;

            p[3] = 0;
            q[3] = 0;
            u[3] = 0;
            finalAdopters[3] = 1;
            initialAdopters[3] = 0;
            initialTime[3] = 0;
            finalTime[3] = 0;

            _setNumber = 0;
            ClearFileData(_path);
        }

        public static void AppendToFile(string stringToAdd, string path = _path)
        {
            File.AppendAllLines(path, new[] { stringToAdd });
        }

        public static void ClearFileData(string path)
        {
            if (File.Exists(path))
                File.Delete(path);
        }

        public static SolPoint[] RK547M(int setNumber)
        {
            SolPoint[] wynik;
            try
            {
                wynik = Ode.RK547M(initialTime[setNumber], initialAdopters[setNumber], (t, x) => (p[setNumber] * (finalAdopters[setNumber] - x)) + ((q[setNumber] * x / (double)finalAdopters[setNumber]) * (finalAdopters[setNumber] - x)) - (u[setNumber] * x)).SolveTo(finalTime[setNumber]).ToArray();
            }
            catch(Exception e)
            {
                Console.WriteLine("Napotkano błąd! Treść: " + e.ToString());
                Console.ReadKey();
                return null;
            }

            return wynik;
        }

        static void Main(string[] args)
        {
            Start();

            string s = "";
            bool loop = true;

            do
            {
                Console.WriteLine("Program implementujący Model Franka Bass'a dotyczący rozprzestrzeniania się \ndanego produktu w czasie, wykorzysujący rozwiązanie ODE.\n\nWartość w nawiasie to aktualna liczba, a w kwadratowym nawiasie zakres.\n\nUWAGA! PROGRAM NIE JEST ODPORNY NA PODAWANIE ZŁYCH ZAKRESÓW/WARTOŚCI!\n");
                Console.WriteLine("\n0 lub Enter - Generuj dane\n1. Współczynnik innowacji - (aktualnie: {0})[0.0 - 1]\n2. Współczynnik imitacji - (aktualnie: {1})[0.0 - 1]\n3. Współczynnik kupców rezygnujących(opcjonalny) - (aktualnie: {2})[0.0 - 1]\n4. Ilość aktualnych kupców - (aktualnie: {3})[1 - {7}]\n5. Ilość ostatecznych kupców - (aktualnie: {4})[1 - {8}]\n6. Czas rozpoczęcia symulacji - (aktualnie: {5})[0 - ?]\n7. Czas zakończenia symulacji - (aktualnie: {6})[0 - ?]\n8. Wybierz zestaw danych (aktualnie: {9})", p[_setNumber].ToString(), q[_setNumber].ToString(), u[_setNumber].ToString(), initialAdopters[_setNumber].ToString(), finalAdopters[_setNumber].ToString(), initialTime[_setNumber].ToString(), finalTime[_setNumber].ToString(), int.MaxValue, int.MaxValue,_setNumber + 1);
                int i;
                int.TryParse(Console.ReadLine(), out i);
                switch (i)
                {
                    default:
                        break;
                    case 0:
                        loop = false;
                        break;
                    case 1:
                        double.TryParse(Console.ReadLine(), out p[_setNumber]);
                        break;
                    case 2:
                        double.TryParse(Console.ReadLine(), out q[_setNumber]);
                        break;
                    case 3:
                        double.TryParse(Console.ReadLine(), out u[_setNumber]);
                        break;
                    case 4:
                        int.TryParse(Console.ReadLine(), out initialAdopters[_setNumber]);
                        break;
                    case 5:
                        int.TryParse(Console.ReadLine(), out finalAdopters[_setNumber]);
                        break;
                    case 6:
                        double.TryParse(Console.ReadLine(), out initialTime[_setNumber]);
                        break;
                    case 7:
                        double.TryParse(Console.ReadLine(), out finalTime[_setNumber]);
                        break;
                    case 8:
                        Console.WriteLine("Wybierz zbiór danych. 1, 2, 3 lub 4, dane zostaną wyświetlone po wybraniu.");
                        int.TryParse(Console.ReadLine(), out _setNumber);
                        _setNumber -= 1;
                        break;
                }
                Console.Clear();
            } while (loop);

            SolPoint[] wynik = RK547M(_setNumber);

            Console.WriteLine("Trwa zapis do pliku...");
            try
            {
                if (wynik != null && wynik.Length > 0)
                {
                    foreach (var item in wynik)
                    {
                        s = "";
                        s += item.X.ToString().Replace("[", string.Empty).Replace("]", string.Empty);
                        s += ";";
                        s += item.T;
                        AppendToFile(s);
                    }
                }
                else
                {
                    Console.WriteLine("Generacja zakończyła się fiaskiem! Złe dane do obliczeń!");
                    return;
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("Napotkano błąd! Treść: " + e.ToString());
                Console.WriteLine("Naciśnij dowolny przycisk, aby zakończyć działanie programu.");
                Console.ReadKey();
                return;
            }
            Console.WriteLine("Pomyślnie zakończono generację pliku danych dla następujących parametrów: \n Współczynnik innowacji = {0}\n Współczynnik imitacji = {1}\nWspółczynnik kupców rezygnujących(opcjonalny) = {2}\n Ilość aktualnych kupców = {3}\n Ilość ostatecznych kupców {4} \n Czas Rozpoczęcia = {5}\n Czas Zakończenia = {6}\nPlik znajduje się w podfolderze projektu pt. 'Generated Files'./nWykres można odtworzyć na podstawie danych z pliku np. w programie GNU Plot.", p[_setNumber], q[_setNumber], u[_setNumber], initialAdopters[_setNumber], finalAdopters[_setNumber], initialTime[_setNumber], finalTime[_setNumber]);
            
            Console.ReadKey();
        }

    }
}
