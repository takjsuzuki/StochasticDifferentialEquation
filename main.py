import Wiener
import OrnsteinUhlenbeck
#import KPZ 

if __name__ == "__main__":

    #Model = Wiener.Wiener()
    Model = OrnsteinUhlenbeck.OrnsteinUhlenbeck()

    Model.Initialize()
    Model.Solve()
    Model.StatisticalQuantities()
    Model.Display()
    Model.Histogram()
    Model.MakeGif()
