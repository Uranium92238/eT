#!/usr/bin/env python3
"""

   eT - a coupled cluster program
   Copyright (C) 2016-2019 the authors of eT
 
   eT is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
 
   eT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.
 
   You should have received a copy of the GNU General Public License
   along with this program. If not, see <https://www.gnu.org/licenses/>.
   This tool is for developers. Run to check any possibile memory leaks:
   it search inside the src directory any call to mem%alloc and mem%dealloc
   and check if all the allocated array have been deallocated with the proper
   dimension. In case of wrong deallocation or not deallcated array, it prints
   the proper warning.
"""

"""

   Check alloc
   Written by Marco Scavino, 2019

   Looks for memory leaks in the source code and identifies possible
   leaks caused by incorrect specification of dimensionalities of arrays
   when allocating/deallocating. This is a script for developers. Note 
   that the script often locates false positives, so each instance must
   be checked manually.

"""

import os
import sys
import re
import argparse
import subprocess

#directories excluded from search
excluded_dir = [
   "tools/linked_list"
]

# get the current repository root directory. If not found, just set empty
repo_dir = subprocess.Popen(['git', 'rev-parse', '--show-toplevel'],
   stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).communicate()[0].rstrip().decode('utf-8')
src_dir = "src"

if sys.version < '3.5':
    print('requires python version >= 3.5')
    sys.exit(1)

def parse_arguments():
   parser = argparse.ArgumentParser()
   parser.add_argument('dir_name', help="directory to check", nargs="?", default="")
   parser.add_argument("-cd", "--curr_dir", help="check in the current directory.", action="store_true")
   args = parser.parse_args()

#  Use the current working directory or the src_dir of the repository
   if args.curr_dir:
      dir_name = os.getcwd()
   else:
      dir_name = os.path.join(repo_dir, src_dir, args.dir_name)

#  chek for the allocation in selected directory
   check_alloc(dir_name)


class bcolors:
   """
   Class for colors in shell
   """
   file_form  = '\033[94m'
   no_dealloc = '\033[92m'
   warning    = '\033[91m'
   reset      = '\033[0m'
   bold       = '\033[1m'
   underline  = '\033[4m'

   def color(self, *formats):
      tmp = ""
      for form in formats:
         tmp += self.__dict__[form]
      return tmp

indent = "   "

class alloc_scan:
   """
   Search all the mem%alloc and mem%dealloc call inside file provide by user.
   """

   def __init__(self):
      self.__alloc_list = []
      self.__alloc_wrong = []
      self.__out = ""
      self.inside_if = 0

   def check_file(self, file_name):
      """
      Check the curent file_name and check the allocation call
      """
      with open(file_name, "r") as file_input:
         self.line = file_input.readline()
         self.line_num = 1
         self.__out = ""
         while self.line:
            self.__check_comment()
            self.line = file_input.readline()
            self.line_num += 1

         if self.__out:
            self.__out = "\n {color1}Wrong deallocation{reset}{out}".format(
               out = self.__out,
               color1 = bcolors.warning,
               reset = bcolors.reset)

         if self.__alloc_list:
            self.__out += "\n {}Not deallocated{}".format(bcolors.warning,bcolors.reset)
            for allocated in self.__alloc_list:
               self.__out += "\n{}{}".format(indent, allocated.get_idx())
      
         if self.__out:
            print("{color}{file_name}{reset}\n{out}\n".format(
               color=bcolors.file_form + bcolors.underline,
               reset=bcolors.reset,
               file_name=file_name,
               out=self.__out))

   def __check_comment(self):
      """
      Check if current line is a comment. If true, skip to the next line
      """
      match = re.search(r'^ *!', self.line)
      if not match:
         match = re.search(r'\bif', self.line)
         if match:
            self.inside_if += 1
         else:
            match = re.search(r'\bend *if\b', self.line)
            self.inside_if -= 1
         self.__check_alloc(self.line)
      return


   def __check_alloc(self, line):
      """
      Check if line contais mem%allocate or mem%deallocate.
      """
      match = re.search(r'call mem%((?:de)?alloc)\((.*?)\)', self.line)
      if match:
         curr_array = array_type(self.line_num, match.group(2).split(","))
         if match.group(1) == "dealloc":
            for idx, allocated in enumerate(self.__alloc_list):
               if allocated.name == curr_array.name:
                  self.__out += allocated.remove_idx(curr_array)
                  if not allocated.indexes:
                     del self.__alloc_list[idx]
                  return
         elif match.group(1) == "alloc":
            for idx, allocated in enumerate(self.__alloc_list):
               if allocated.name == curr_array.name:
                  allocated.append_idx(curr_array.indexes[0])
                  return
            self.__alloc_list.append(curr_array)
      return

class array_type:
   """
   Class to save the array name, the line of allocation and the indexes
   used for allocated the array.
   """

   def __init__(self, line:int, indexes:str):
      self.name = indexes[0]
      self.line = [line]
      tmp = []
      for idx in indexes[1:]:
         tmp.append("".join(idx.split()))
      self.indexes = [tmp]

   def __eq__(self, array):
      return self.name == array.name and self.indexes == array.indexes
   
   def get_idx(self):
      tmp = ""
      for (idx, line) in zip(self.indexes, self.line):
         tmp += "{color1}{array}({idx}){reset} in line {color2}{line}{reset}".format(
            color1=bcolors.bold,
            color2=bcolors.warning,
            reset=bcolors.reset,
            array=self.name,
            idx=" ".join(idx),
            line=line)
      return tmp

   def append_idx(self, idx:list):
      self.indexes.append(idx)

   def remove_idx(self, curr_idx) -> str:
      try:
         self.indexes.remove(curr_idx.indexes[0])
         return ""
      except ValueError:
         return "\n{0}{1} deallocated as\n{0}{2}\n".format(indent, self.get_idx(), curr_idx.get_idx())

def check_alloc(dir_name:str, ext:str=".F90"):
   file_list = []
   for (dir_path, dir_names, file_names) in os.walk(dir_name):
         for file_name in file_names:
            curr_ext = os.path.splitext(file_name)[1]
            exclude = False
            for dir_ in excluded_dir:
               if dir_ in dir_path:
                  exclude = True
                  break
            if curr_ext == ext and not exclude:
               alloc_dir = alloc_scan()

               alloc_dir.check_file(os.path.join(dir_path, file_name))

if __name__ == "__main__":
   parse_arguments()