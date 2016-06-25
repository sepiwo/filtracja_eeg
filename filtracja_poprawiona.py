import matplotlib.pyplot as plt
import pandas
from numpy import *
from scipy.signal import iirdesign, freqz, iirfilter, lfilter, group_delay


#data_file="D:/STUDIA/Projekty/EEG/Projekt EEG/skrypty/csv/BCICIV_2a_gdf/A01T_signal.csv"
#info_file="D:/STUDIA/Projekty/EEG/Projekt EEG/skrypty/csv/BCICIV_2a_gdf/A01T_info.csv"
data_file="/home/sebastian/studia/projekt/eeg/dane/A01T_signal.csv"

info_file="/home/sebastian/studia/projekt/eeg/dane/A01T_info.csv"


class filter_parameters: #parametry filtru

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
        #self.num, self.denom = iirdesign(wp=[self.filter_parameters.lowcut, self.filter_parameters.highcut], ws=[0.01,0.35], gpass=0.05, gstop=5, analog=False, ftype='cheby1', output='ba')
        #self.num, self.denom = iirdesign(wp=[self.filter_parameters.lowcut, self.filter_parameters.highcut], ws=[0.005,0.35], gpass=1.7, gstop=20, analog=False, ftype='cheby2', output='ba')
        #self.num, self.denom = iirdesign(wp=[self.filter_parameters.lowcut, self.filter_parameters.highcut], ws=[0.005,0.35], gpass=0.4, gstop=28, analog=False, ftype='ellip', output='ba')
        self.num, self.denom = iirfilter(5, [self.filter_parameters.lowcut, self.filter_parameters.highcut],btype='bandpass', analog=False, ftype='bessel', output='ba')


def pair(data, labels=None):
    """ Generate something similar to R `pair` """

    nVariables = data.shape[1]
    if labels is None:
        labels = ['var%d'%i for i in range(nVariables)]
    fig = plt.figure()
    for i in range(nVariables):
        for j in range(nVariables):
            nSub = i * nVariables + j + 1
            ax = fig.add_subplot(nVariables, nVariables, nSub)
            if i == j:
                ax.hist(data[:, i])
                ax.set_title(labels[i])
            else:
                ax.plot(data[:, i], data[:, j], '.k')

    return fig



def replace_nan(dane):
        dane[isnan(dane)]=0
        return dane

def signal_energy(przefiltrowany,window_size):
    przefiltrowany=square(przefiltrowany)
    energy= przefiltrowany
    
    for i in range(0,energy.shape[0]-1):
        for j in range(0,(energy.shape[1]-window_size-1)):
            energy[i,j+window_size-1]=mean(przefiltrowany[i,j:(j+window_size)])
    return energy


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
        self.energy=self.raw_signal

    def description(self): #wyswietlanie informacji o sygnalach
        print('Signal inormations:\n channels: ' ,self.raw_signal.shape[1] ,'\n signal length:' ,self.raw_signal.shape[0], '\n classes: ', pandas.unique(self.event) )

    def plot_raw(self): #wyswietlenie wykresu dla danych surowych dla parametrow ustawionych przy pomocy .choosen_channels oraz .start_time i .end_time

        plt.plot(self.raw_signal.ix[:,self.choosen_channels],'b-', label = 'Sygnal nieprzefiltroway')
        plt.title('Raw EEG signal')
        plt.xlabel('Numer probki')
        plt.ylabel("Amplituda sygnalu EEG")


    def perform_filtration(self):#filtracja sygnalu!
        processed=empty([self.raw_signal.shape[0],self.raw_signal.shape[1]])
        for i in range(0,(self.raw_signal.shape[1]-1)):
            processed[:,i]=lfilter(self.freq_filter.num, self.freq_filter.denom, self.raw_signal.ix[:,1])
        return processed

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
        znie = plt.plot((self.filter_parameters.sample_rate*0.5/pi)*w,angles, 'g-', label = 'Char. fazowa')
        plt.ylabel('Faza')

        plt.grid()
        tekst = powa + znie
        wybierz = [l.get_label() for l in tekst]

        plt.legend(tekst, wybierz, loc='best')

        #####   opoznienie grupowe    ##########
        plt.subplot(2,1,2)

        w2, gd = group_delay((self.freq_filter.num, self.freq_filter.denom))

        plt.plot((self.filter_parameters.sample_rate*0.5/pi)*w2, gd)
        plt.grid()
        plt.xlabel('Czestotliwosc [Hz]')
        plt.ylabel('Opoznienie grupowe [probki]')
        plt.title('Opoznienie grupowe filtru')

        plt.show()





#pobieranie danych
data=eeg_data(data_file,info_file)

#wyswietlanie informacji o sygnalach
data.description()

#wyswietlanie sygnalow wszystkich
plot=data.plot_raw()
plt.show()


                                  ########### RYSOWANIE  CHARAKTERYSTYK FILTRU ################
data.plot_filter_characteristics()

#wyswietlanie konkretnego sygnalu
data.choosen_channels=1
wykres = data.plot_raw()
plt.title("Sygnal surowy, kanal 1-szy ")
plt.show()


                                ################   FILTROWANIE     ################
przefiltrowany=data.perform_filtration()

plt.plot(przefiltrowany[:,1], 'g-', label = 'Sygnal przefiltrowany')
plt.title('Sygnal przefiltrowany')
plt.xlabel('Numer probki')
plt.ylabel('Amplituda sygnalu EEG')
plt.show()


####Energia

energy=signal_energy(przefiltrowany,50)
plt.plot(energy[:,1])


plt.show()

#przefiltrowany.shape

#przefiltrowany
#data.raw_signal[:,1]

# obiekty opisane przez 5 cech (energie dla 5-ciu kanalow)
dane_5cech = zeros((len(energy[:, 3]), 5))
dane_5cech[:, 0] = energy[:, 0]  # pamietajac ze Python liczy od zera (wiec nr kanalu to + 1)
dane_5cech[:, 1] = energy[:, 7]
dane_5cech[:, 2] = energy[:, 9]
dane_5cech[:, 3] = energy[:, 11]
dane_5cech[:, 4] = energy[:, 19]
etykiety = ["kanal 1", "kanal 8", "kanal 10", "kanal 12", "kanal 20"]
pair(dane_5cech, etykiety)
plt.show()

# przefiltrowany.shape

# przefiltrowany
# data.raw_signal[:,1]

