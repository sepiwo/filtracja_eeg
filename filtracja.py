
import matplotlib.pyplot as plt
import pandas
from numpy import *
from scipy.signal import iirdesign, freqz, iirfilter, lfilter
from scipy.signal import freqz, group_delay
from matplotlib.pyplot import *





#data_file=pandas.read_csv(raw_input('path to signal CSV file: '))
#info_file=pandas.read_csv(raw_input('path to info CSV file: '))

data_file="/home/sebastian/studia/projekt/eeg/dane/A01T_signal.csv"

info_file="/home/sebastian/studia/projekt/eeg/dane/A01T_info.csv"

class eeg_data :

    def __init__(self,data_file,info_file): #tworzenie obiektu przechowujacego sygnaly oraz informacje o nim
        self.raw_signal=pandas.read_csv(data_file)  #surowy sygnal
        self.event=pandas.read_csv(info_file).ix[:,3] #wektor zawierajacy informacje o zdarzeniach w czasie trwania sygnalu
        self.start_time=pandas.read_csv(info_file).ix[:,1]  #wektor zawierajacy czas poczatkowy
        self.end_time=pandas.read_csv(info_file).ix[:,2]  #wektor zawierajacy czas koncowy
        self.time=pandas.read_csv(info_file).ix[:,0]    #
        self.choosen_channels=range(0,self.raw_signal.shape[1]-1) #wybrane sygnaly do sporzadzenia wykresu domyslnie wszystkie!
        self.choosen_time_start=0
        self.choosen_time_stop=1
        self.choosen_events=[pandas.unique(self.event)] #wybrane zdarzenia domyslnie wszystkie
        self.filtered=self.raw_signal    #sygnal przefitrowany
        self.krok_czasu=self.time[1]-self.time[0]



    def description(self): #wyswietlanie informacji o sygnalach
        print('Signal inormations:\n channels: ' ,self.raw_signal.shape[1] ,'\n signal length:' ,self.raw_signal.shape[0], '\n classes: ', pandas.unique(self.event) )
    def plot_raw(self): #wyswietlenie wykresu dla danych surowych dla parametrow ustawionych przy pomocy .choosen_channels oraz .start_time i .end_time
        plt.plot(self.raw_signal.ix[self.choosen_channels,:])
    def filter(self,):
        pass



data=eeg_data(data_file,info_file)

data.description()


data.plot_raw()

data.filter()

#8,10,12
#23,24,25 <-EOG

#POS pozycja (w ktorej probce zapiasna  chwila czasowa)
#Event_typ co sie wydarzylo -> odczytujemy z artykulu w tabeli
#EVENT.DUR jak cos dlugo trwalo pos+dur czas do kiedy
#wyobrazanie sobie od 3 do 6 sek
#dane do analizy h.EVENT.POS(n)+3*probkowanine: h.EVENT.POS(n)+6*probkowanine
#wydobyc 288 wyobrazen sugnalow o dlugosci 3 sekund (dzielimy nasz sygnal na klasy)
#
#Wczytywanie i dzielenie sygnalow <SEBASTIAN>
#fILTRACJA PRZESTRZENNA
#



def char_plot(num, denom, sample_rate):

    w, h = freqz(num, denom)
    figure(1)
    subplot(2,1,1)
    hold(True)
    pow = plot((sample_rate*0.5/pi)*w, abs(h),'b-', label = 'Char. amplitudowa')
    title('Charakterystyki filtru')
    xlabel('Czestotliwosc [Hz]')
    ylabel('Amplituda')


    twinx(ax=None)
    angles = unwrap(angle(h))
    aznie = plot((sample_rate*0.5/pi)*w, angles, 'g-', label = 'Char. fazowa')
    ylabel('Faza')

    grid()
    tekst = pow + aznie
    wybierz = [l.get_label() for l in tekst]

    legend(tekst, wybierz, loc='best')
########################################################################################################################
    subplot(2,1,2)

    w2, gd = group_delay((num, denom))

    plot((sample_rate*0.5/pi)*w2, gd)
    grid()
    xlabel('Czestotliwosc [Hz]')
    ylabel('Opoznienie grupowe [probki]')
    title('Opoznienie grupowe filtru')

    show()









close("all")

sample_rate = 250 #Hz
f1 = 8 #8Hz
f2 = 30 #30Hz

nyq = 0.5*250       #czestotliwosc Nyquista
lowcut = f1/nyq     #skalowanie do czestoliwosci Nyquista
highcut = f2/nyq



                                                #################   BUTTER ################

mloss = 5
matten = 20

num, denom = iirdesign(wp=[lowcut, highcut], ws=[0.005,0.35], gpass=mloss, gstop=matten, analog=False, ftype='butter', output='ba')



                                                #################   CHEBY1  ################

mloss2 = 0.05
matten2 = 5

num2, denom2 = iirdesign(wp=[lowcut, highcut], ws=[0.01,0.35], gpass=mloss2, gstop=matten2, analog=False, ftype='cheby1', output='ba')


                                                #################   CHEBY2  ################

mloss3 = 1.7            #tu pozmienialam wartosci, zeby zmniejszyc wielkosci tetnien w pasmie zaporowym
matten3 = 20

num3, denom3 = iirdesign(wp=[lowcut, highcut], ws=[0.005,0.35], gpass=mloss3, gstop=matten3, analog=False, ftype='cheby2', output='ba')



                                                #################   ELIPTYCZNY  ################

mloss4 = 0.4            #tu tez pozmienialam wartosci
matten4 = 28

num4, denom4 = iirdesign(wp=[lowcut, highcut], ws=[0.005,0.35], gpass=mloss4, gstop=matten4, analog=False, ftype='ellip', output='ba')



                                                #################   BESSEL  ################

#iirfilter, bo iirdesign nie daje rady z Besselem...

num5, denom5 = iirfilter(5, [lowcut, highcut],btype='bandpass', analog=False, ftype='bessel', output='ba')


                                            ########### RYSOWANIE   ################

char_plot(num5, denom5, sample_rate)

                                ################   FILTROWANIE     ################
przefiltrowany = lfilter(num, denom, data.raw_signal.ix[:,0])
