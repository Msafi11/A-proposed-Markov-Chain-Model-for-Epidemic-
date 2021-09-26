import numpy as np
import random as rm
import quantecon as qe
import tkinter
from tkinter import *
from tkinter import simpledialog
from tkinter import messagebox
  

root = tkinter.Tk()
root.geometry("300x300")
root.withdraw()
###################################################################################
def get_trans_matrix():   
    #Function to get the probabilities and make the transistion matrix
    p = np.zeros((3,3),dtype=float)  # Initialize the transistion matrix
    for i in range(3):
        pi = p[i]
        for j in range(3):
            
            pi[j] = simpledialog.askfloat(title="MarkovChain Prediction",prompt="Enter P" + str(i+1) + str(j+1) + " : ",minvalue=0.0, maxvalue=1.0,parent=root)
           
            
        if((np.sum(p[i])) != 1.0): 
            while(np.sum(p[i]) != 1.0):
                messagebox.showinfo("MarkovChain Prediction","Wrong Sum of all Pi must be equal  1.!.!.Re Enter it please:")
                
                for j in range(3):
                    pi[j] = simpledialog.askfloat(title="MarkovChain Prediction",prompt="Enter P" + str(i+1) + str(j+1) + " : ",minvalue=0.0, maxvalue=1.0)
                                   
           
    return p
    

###############################################################################################
def check_MC_reducible(transitionmatrix): #check this markov chain is irreducible or reducible
    mc = qe.MarkovChain(transitionmatrix, ('infected','recovered', 'dead'))
    if(mc.is_irreducible):
        messagebox.showinfo("MarkovChain Prediction","This MC is irreducible")
    else:
        messagebox.showinfo("MarkovChain Prediction","This MC is reducible")  
############################################################################################
def check_periodicity(transitionmatrix): # check the periodicity of the markov chain and the states
    mc = qe.MarkovChain(transitionmatrix)
    if(mc.is_aperiodic):
        messagebox.showinfo("MarkovChain Prediction","This MC is aperiodic and all states are aperiodic ")
        
    else:
        period = mc.period
        messagebox.showinfo("MarkovChain Prediction","This MC is periodic and all states are aperiodic with period = " + str(period)) 
#############################################################
    
def get_reachable_states(transitionmatrix): # Getting The reachable states and the not reachable
    for i in range(3): 
        pi = transitionmatrix[i] 
        for j in range(3):
            if(pi[j] > 0.0): # A state j is said to be [reachable] from a state i if Pij > 0 for some n = 1, 2, ...
                messagebox.showinfo("MarkovChain Prediction","State " + str(j) + " is reachable from state "+ str(i))
            else:
                messagebox.showinfo("MarkovChain Prediction","State " + str(j) + " is not reachable from state "+ str(i)) 
#####################################################################################
def get_absorbing_states(transitionmatrix): # Getting The absorbing states and the not absorbing
    for i in range(3):
        pi = transitionmatrix[i]
        if(pi[i] == 1.0): #A state i is said to be absorbing if it forms a single element closed set. If i is an absorbing state we have Pij =1
            messagebox.showinfo("MarkovChain Prediction","State " + str(i) + " is an absorbing state")
        else:
            messagebox.showinfo("MarkovChain Prediction","State " + str(i) + " is not an absorbing state")
###########################################################################
def get_recurrent_transient_states(transitionmatrix): # Getting The recurrent states and the transient
    for i in range(3):
        pi = transitionmatrix[i]
        if(pi[i] == 1.0): # A state i is said to be recurrent if Pi = 1
            messagebox.showinfo("MarkovChain Prediction","State " + str(i) + " is recurrent state")
        elif(pi[i] < 1.0): # and if Pi < 1, The State is transient state
            messagebox.showinfo("MarkovChain Prediction","State " + str(i) + " is transient state")
#############################################################################

def prediction(n,bi0,transitionmatrix,states): # Function to predict the next state
    # n : the period we want do predect the next state after it(n days, n months, n years...etc)
    # bi0 : The intial Bi 
    # states : The states names
    transitionmatrix_n = transitionmatrix
    if(n == 1): # if the period is 1 so we already have the transitionmatrix of P 
        transitionmatrix_n = transitionmatrix
    else: # if the Period more than 1 so we need to calculate the P power n
        for i in range(n):
            if(i == 0):
                transitionmatrix_n = transitionmatrix
            else:
                transitionmatrix_n = np.matmul(transitionmatrix_n, transitionmatrix_n)        
    messagebox.showinfo("MarkovChain Prediction","P after "+str(n)+ " days is:\n "+str(transitionmatrix_n))
    bi_n = np.matmul(bi0, transitionmatrix_n) # calculate Bi power N
    messagebox.showinfo("MarkovChain Prediction","bi after "+str(n)+ " days is: " + str(bi_n))
    messagebox.showinfo("MarkovChain Prediction","The next state of the virus after " + str(n) + " days will be: ")
    for i in range(3):
        messagebox.showinfo("MarkovChain Prediction",str(i+1)+"- " + states[i] + " with probability " + str(bi_n[i]))
    messagebox.showinfo("MarkovChain Prediction","So with this probabilities we can predict that the number of " + states[np.argmax(bi_n)] + " people will be increased in the next " + str(n) + " Days") #Predict the next state 
#######################################################################################

states = ["infected","recovered","dead"] # states names
bi0 = np.zeros(1) 
bi0 = [1.0,0.0,0.0] # Initial Bi        
p  = get_trans_matrix() #Getting the transistion Matrix
messagebox.showinfo("MarkovChain Prediction","The Transition matrix is:\n"+str(p))
check_MC_reducible(p)
check_periodicity(p)
get_absorbing_states(p)
get_reachable_states(p)
get_recurrent_transient_states(p)
n = simpledialog.askinteger(title="MarkovChain Prediction",prompt="Enter The period of time you want to predict the state after it")
prediction(n,bi0,p,states)
#######################################################################################
