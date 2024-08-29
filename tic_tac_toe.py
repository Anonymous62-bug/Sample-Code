#!/usr/bin/python3

### This is a program to play tic-tac-toe
### Requirements: Python3 and NumPy,Pandas packages
### To run:python3 filename.py

import numpy as np
import pandas as pd

class tictactoe:
    player = "A"
    board_matrix = np.zeros(9).reshape(3,3)
    game_complete = False
    player_values = {'A': 1, 'B': 4}

    #initialization
    def __init__(self):
        print("Let's play TicTacToe!\n")
        print("Player A is X's\nPlayer B is O's")

    #Gets coordinates from player
    #Returns: (str) User input
    def get_coordinates(self):
        user_input = input("Location in the tictactoe box:").rstrip("\n")
        return(user_input)

    #Checks if user input is valid
    #Return: (int,int) x and y coordinates
    def check_input(self,user_input):
        try:
            x_coord,y_coord = user_input.split(",")
            x_coord,y_coord = int(x_coord),int(y_coord)
            for coord in [x_coord,y_coord]:
                if coord not in [0,1,2]:
                    print("Box location has to be between 0,1 or 2, try again!")
                    return(self.check_input(self.get_coordinates()))
                else:
                    return(x_coord,y_coord)
        except ValueError as err:
            print("ValueError {}".format(err))
            return(self.check_input(self.get_coordinates()))

    #Outputs the game board
    #Return: NULL
    def print_box(self):
        print_matrix = np.array(['']*9).reshape(3,3)
        print_matrix[self.board_matrix == 0] = "-"
        print_matrix[self.board_matrix == 1] = "X"
        print_matrix[self.board_matrix == 4] = "O"
        board = """
           0 1 2
        0  {}|{}|{}
           ------
        1  {}|{}|{}
           ------
        2  {}|{}|{}
        """.format(print_matrix[0,0],print_matrix[0,1],print_matrix[0,2],
        print_matrix[1,0],print_matrix[1,1],print_matrix[1,2],
        print_matrix[2,0],print_matrix[2,1],print_matrix[2,2])
        print(board)

    #Checks if value can be placed at coordinate
    #Return: (bool)
    def check_valid_insert(self,x_coord,y_coord):
        if self.board_matrix[x_coord,y_coord] == 0:
            return(True)
        else:
            print("Already inserted, pick another spot!")
            return(False)
    #Inserts the value at coordinate
    #Return: NULL
    def insert_input(self,x_coord,y_coord):
        self.board_matrix[x_coord,y_coord] = self.player_values[self.player]
        return()

    #Checking if game is over (win or tie)
    #Returns: NULL
    def check_game_status(self):
            board_df = pd.DataFrame(self.board_matrix)
            sum_diag1,sum_diag2 = sum(np.diag(self.board_matrix)),sum(np.diag(np.fliplr(self.board_matrix)))
            for key in self.player_values.keys():
                check_status = list()
                #checking columns
                check_status.append(True in list(board_df.prod(axis=1) == self.player_values[key]**3))
                #checking rows
                check_status.append(True in list(board_df.prod(axis=0) == self.player_values[key]**3))
                #checking diagonals
                check_status.append(sum_diag1 == self.player_values[key]*3 or sum_diag2 == self.player_values[key]*3)
                if True in check_status:
                    print("Player {} wins!".format(key))
                    self.game_complete = True
                    self.print_box()
                    return()
            if 0 not in self.board_matrix:
                print("Game Over, No one wins :(")
                self.game_complete = True
                self.print_box()
                return()
            else:
                return()
#Main function
def main():
    play = tictactoe()
    while not play.game_complete:
        play.print_box()
        print("Player {}'s turn".format(play.player))
        x,y = play.check_input(play.get_coordinates())
        if play.check_valid_insert(x,y):
            play.insert_input(x,y)
            play.check_game_status()
            if play.player == "A":
                play.player = "B"
            else:
                play.player = "A"

if __name__ == "__main__":
    main()
