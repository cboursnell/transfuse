require 'open3'

module Transfuse

  class Status
    def success?
      return true
    end
  end

  class Cmd

    attr_accessor :cmd, :stdout, :stderr, :status

    def initialize cmd
      @cmd = cmd
    end

    def run file=nil
      unless file.nil?
        if File.exist?(file)
          @stdout = ""
          @stderr = ""
          @status = Status.new
          return true
        end
      end
      @stdout, @stderr, @status = Open3.capture3 @cmd
      return false
    end

    def to_s
      @cmd
    end

  end

end
