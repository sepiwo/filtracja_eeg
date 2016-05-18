

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



class filter_parameters: #parametry filtru Czestotliwosciowego
    def __init__(self,sample_rate=250,f1=8,f2=30,nyq=0.5*250):
        self.sample_rate = sample_rate #Hz #czestotliwosc probkowania
        self.f1 = f1 #8Hz    czestotliwosc dolna
        self.f2 = f2 #30Hz czestotliwosc gorna
        self.nyq = nyq       #czestotliwosc Nyquista
        self.lowcut = self.f1/self.nyq     #skalowanie do czestoliwosci Nyquista
        self.highcut = self.f2/self.nyq
        self.mloss = 5
        self.matten = 20
        self.type="bessel"



class freq_filter:
    def __init__(self):
        self.filter_parameters=filter_parameters()
        #self.num, self.denom = iirdesign(wp=[self.filter_parameters.lowcut, self.filter_parameters.highcut], ws=[0.005,0.35], gpass=self.filter_parameters.mloss, gstop=self.filter_parameters.matten, analog=False, ftype='butter', output='ba')
        self.num, self.denom = iirfilter(5, [self.filter_parameters.lowcut, self.filter_parameters.highcut],btype='bandpass', analog=False, ftype='bessel', output='ba')



def replace_nan(dane):
        dane[isnan(dane)]=0
        return dane



class eeg_data :

    def __init__(self,data_file,info_file): #tworzenie obiektu przechowujacego sygnaly oraz informacje o nim
        self.raw_signal=replace_nan(pandas.read_csv(data_file))  #surowy sygnal
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
        self.filter_parameters=filter_parameters()
        self.freq_filter=freq_filter()

    def description(self): #wyswietlanie informacji o sygnalach
        print('Signal inormations:\n channels: ' ,self.raw_signal.shape[1] ,'\n signal length:' ,self.raw_signal.shape[0], '\n classes: ', pandas.unique(self.event) )
    def plot_raw(self): #wyswietlenie wykresu dla danych surowych dla parametrow ustawionych przy pomocy .choosen_channels oraz .start_time i .end_time
        plt.plot(self.raw_signal.ix[:,self.choosen_channels])
        plt.title('Raw EEG signal')
        plt.xlabel('Time [ms]')
    def perform_filtration(self):#filtracja sygnalu!
        return lfilter(self.freq_filter.num, self.freq_filter.denom, data.raw_signal.ix[:,4])

    def plot_filter_characteristics(self):
        w, h = freqz(self.freq_filter.num, self.freq_filter.denom)
        plt.figure(1)
        plt.subplot(2,1,1)
        plt.hold(True)
        powa = plt.plot((self.filter_parameters.sample_rate*0.5/pi)*w, abs(h),'b-', label = 'Char. amplitudowa')
        plt.title('Charakterystyki filtru')
        plt.xlabel('Czestotliwosc [Hz]')
        plt.ylabel('Amplituda')


        plt.twinx(ax=None)
        angles = unwrap(angle(h))
        plt.znie = plot((self.filter_parameters.sample_rate*0.5/pi)*w,angles, 'g-', label = 'Char. fazowa')
        plt.ylabel('Faza')

        plt.grid()
        tekst = powa + znie
        wybierz = [l.get_label() for l in tekst]

        plt.legend(tekst, wybierz, loc='best')
    ########################################################################################################################
        plt.subplot(2,1,2)

        w2, gd = group_delay((num, denom))

        plt.plot((sample_rate*0.5/pi)*w2, gd)
        plt.grid()
        plt.xlabel('Czestotliwosc [Hz]')
        plt.ylabel('Opoznienie grupowe [probki]')
        plt.title('Opoznienie grupowe filtru')

        plt.show()






data=eeg_data(data_file,info_file)

data.description()
unique(isfinite(data.raw_signal))



plot=data.plot_raw()
plt.show()

data.plot_filter_characteristics()


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
#przefiltrowany = lfilter(num, denom, data.raw_signal.ix[:,0])
przefiltrowany=data.perform_filtration()
plt.plot(przefiltrowany)
plt.show()


subplot(1,1,1)
