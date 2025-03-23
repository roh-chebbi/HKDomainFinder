#!/usr/bin/env python

"""usage: ./rchebbi3_question1.py <character> <height>"""

import sys

def printTriangle(char_type:str, calculated_args: tuple) -> None:
    """Prints triangle based on character and height
    
    Args: 
        char_type: user input
        calculated_args: tuple of calculated new height and if the height is odd or even
    
    Returns:
        None
    """
    
    new_height = calculated_args[0]
    flag       = calculated_args[1]

    for i in range(1,new_height,1):
        print(i*char_type)

    if flag == "even":
        print(new_height*char_type)

    for i in range(new_height,0,-1):
        print(i*char_type)

def calculateHeight(height: int) -> tuple:
    """Calculates if height is odd or even
    
    Args: 
        height: user input
    
    Returns: 
        New height and flag for odd/even
    """
    
    if height%2 == 0:
        new_height = int(height/2)
        flag = "even"
    else:
        new_height = int(height/2) + 1
        flag = "odd"
    
    return(new_height,flag)


calculated_args=calculateHeight(int(sys.argv[2]))
printTriangle(sys.argv[1],calculated_args)